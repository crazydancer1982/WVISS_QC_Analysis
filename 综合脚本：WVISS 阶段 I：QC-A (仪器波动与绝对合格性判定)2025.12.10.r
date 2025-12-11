# ======================================================================
# 综合脚本：WVISS 阶段 I：QC-A (仪器波动与绝对合格性判定)
# ----------------------------------------------------------------------
# 包含：数据加载、预处理、Cycle-Level 统计、Z-Composite 计算和 CL/SL 判定
# ----------------------------------------------------------------------

# 0. 环境设置与库加载
# 确保已安装：dplyr, tidyr, zoo, RcppRoll
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
if (!requireNamespace("RcppRoll", quietly = TRUE)) install.packages("RcppRoll")

library(dplyr)
library(tidyr)
library(zoo)
library(RcppRoll)

# 设置有效数字位数
options("digits" = 9) 

# --- 全局 QC-A 参数设定 (采用最终折衷方案的推荐值) ---
N_WINDOW <- 270     # 窗口长度 N* (Cycles)
K_CL_SIGMA <- 3     # 动态 CL 阈值 K (K*MAD)
K_STEP <- 5         # 步长 k* (Cycles)
K_SL_SIGMA <- 5     # 静态 SL 阈值 K (K*MAD)

# 0.1 定义稳健性 MAD 函数
robust_mad <- function(x) {
  # 确保在计算 median 和 mad 时处理 NA
  mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
}


# 0.2 定义稀疏滚动 MAD 函数 (用于高效计算 CL)
calculate_sparse_rolling_mad <- function(data_vector, window_n, k_step) {
  L <- length(data_vector)
  sparse_mad <- rep(NA_real_, L)
  
  # 确定计算的索引：从 N_WINDOW 开始，每隔 k_step 计算一次
  start_idx <- N_WINDOW

  indices_to_calc <- seq(start_idx, L, by = k_step)

  for (i in indices_to_calc) {
    window_data <- data_vector[(i - window_n + 1):i]
    mad_val <- robust_mad(window_data)
    sparse_mad[i] <- mad_val
  }

  # 【增强逻辑 1: 处理序列开始 (1 到 window_n - 1)】
  if (L >= window_n && window_n > 1) {
    first_mad_val <- sparse_mad[window_n]
    # 采用第一个完整窗口的值进行回填
    if (!is.na(first_mad_val)) {
      sparse_mad[1:(window_n - 1)] <- first_mad_val
    }
  }

  # 使用 Last Observation Carried Forward (na.locf) 填充 NA
  sparse_mad <- zoo::na.locf(sparse_mad, na.rm = FALSE)
  return(sparse_mad)
}

# ======================================================================
# 1. 数据加载与预处理
# ======================================================================

# 设置工作路径
wd.synthesis='/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'
setwd(wd.synthesis)

#读如WVISS data
WVISS.all=read.csv('WVISS.all.csv',header=T, stringsAsFactors = FALSE)

message("WVISS.all.csv 数据加载完成。")

# --- 1.1 核心变量计算 ---
RD_smow=0.00015576 # VSMOW of D
RO_smow=0.0020052  # VSMOW of 18O

# 还原 DOY_mean 为日期
WVISS.all$date=as.POSIXct(86400* WVISS.all$DOY_mean, origin = "2011-01-01")

# 计算观测值 dO_obs 和 dD_obs
WVISS.all$dO_obs=((WVISS.all$O18ppm_mean/WVISS.all$O16ppm_mean)/RO_smow-1)*1000
WVISS.all$dD_obs=((WVISS.all$HODppm_mean/(WVISS.all$O16ppm_mean*2))/RD_smow-1)*1000

# 计算原始偏差 (raw bias)
WVISS.all$dO_raw_bias=WVISS.all$dO_obs-WVISS.all$o_wahaha
WVISS.all$dD_raw_bias=WVISS.all$dD_obs-WVISS.all$d_wahaha

# 创建唯一观测值标识
WVISS.all$group.gradient=paste(WVISS.all$group.id,WVISS.all$Gradient,sep='_')


# --- 1.2 按 Cycle 时间（group_DOY_mean）进行排序 ---
# 计算每个 Cycle 的平均 DOY
group_DOY_mean=WVISS.all %>%
  group_by(group.id) %>%
  summarise(group_DOY_mean=mean(DOY_mean, na.rm = TRUE), .groups = 'drop')

# 合并平均 DOY 并按时间排序
WVISS.all=merge(WVISS.all,group_DOY_mean,by='group.id')
WVISS.all=WVISS.all[order(WVISS.all$group_DOY_mean),]

# --- 1.3 移除冗余或不需要的列 ---
WVISS.all=WVISS.all[,-which(names(WVISS.all) %in% c('SD_outlier','keep.status','month', 'X.1','X','DOY_SD'))]

# --- 1.4 Cycle 级别中心化偏差计算 (dD_center, dO_center) ---
WVISS.all = WVISS.all %>%
    group_by(group.id) %>%
    mutate(
        # dD_center = 原始偏差 - group内平均原始偏差
        dD_center = dD_raw_bias-mean(dD_raw_bias, na.rm = TRUE),
        dO_center = dO_raw_bias-mean(dO_raw_bias, na.rm = TRUE)
    ) %>%
    ungroup()

# --- 1.5 水汽浓度 (H2O_mean) 全局 Z-Score 标准化 (H_Zscore) ---
H_mean_global <- mean(WVISS.all$H2O_mean, na.rm = TRUE)
H_sd_global <- sd(WVISS.all$H2O_mean, na.rm = TRUE)

WVISS.all <- WVISS.all %>%
    mutate(
        # H_Zscore 是基于全局均值和标准差的标准化变量
        H_Zscore = (H2O_mean - H_mean_global) / H_sd_global
    )
message("数据预处理和中心化计算完成。")


# ======================================================================
# 2. QC-A 核心指标计算：Cycle-Level 统计与 Z-Composite
# ======================================================================

# --- 2.1 Cycle-Level 变异性指标计算 (QC-A 基础数据) ---
df_cycle_metrics <- WVISS.all %>%
  # 确保数据仍按 group_DOY_mean 排列
  arrange(group_DOY_mean) %>%
  group_by(group.id, group_DOY_mean, year) %>%
  summarise(
    # 变异性指标：Cycle 内中心化偏差的标准差 (SD)
    dD_center_SD = sd(dD_center, na.rm = TRUE),
    dO_center_SD = sd(dO_center, na.rm = TRUE),
    # 绝对合格指标：H2O Z-score 的均值（绝对值）
    H_Zscore_Mean_abs = abs(mean(H_Zscore, na.rm = TRUE)),
    .groups = 'drop'
  ) %>%
  # 添加 Cycle 序列号 (QC 判定基于此序列)
  mutate(Cycle_Seq = 1:n())

# --- 2.2 计算全局稳健统计量 (Median 和 MAD) ---
global_stats <- df_cycle_metrics %>%
  summarise(
    D_SD_Median = median(dD_center_SD, na.rm = TRUE),
    O18_SD_Median = median(dO_center_SD, na.rm = TRUE),
    H_Z_Median = median(H_Zscore_Mean_abs, na.rm = TRUE),
    
    D_SD_MAD = robust_mad(dD_center_SD),
    O18_SD_MAD = robust_mad(dO_center_SD),
    H_Z_MAD = robust_mad(H_Zscore_Mean_abs)
  )

# --- 2.3 标准化：计算 Z-Scores 和 Z-Composite (综合 Z-Score) ---
df_qc_a <- df_cycle_metrics %>%
  mutate(
    # Cycle SD Z-Score
    Z_SD_D = (dD_center_SD - global_stats$D_SD_Median) / global_stats$D_SD_MAD,
    Z_SD_O18 = (dO_center_SD - global_stats$O18_SD_Median) / global_stats$O18_SD_MAD,
    # H2O 绝对 Z-Score
    #Z_H2O_abs = (H_Zscore_Mean_abs - global_stats$H_Z_Median) / global_stats$H_Z_MAD,
    
    # Z-Composite = max(Z_SD, Z_H2O_abs) (取最大值作为最差情况)
    Z_Composite_D = Z_SD_D
    Z_Composite_O18 = Z_SD_O18
  )

head(as.data.frame(df_qc_a))
# ======================================================================
# 3. 阈值计算与 Outlier 判定
# ======================================================================

# --- 3.1 动态控制限 (CL) 计算 ---
# CL = Rolling Median + K_CL_SIGMA * Rolling MAD
df_qc_a <- df_qc_a %>%
  mutate(
    # 滚动中位数 (Rolling Median) - 用于趋势拟合
    Rolling_Median_D = roll_median(Z_Composite_D, n = N_WINDOW, align = "right", fill = NA, na.rm = TRUE),
    Rolling_Median_O18 = roll_median(Z_Composite_O18, n = N_WINDOW, align = "right", fill = NA, na.rm = TRUE)
  ) %>%
  mutate(
    # 稀疏滚动 MAD (Rolling MAD) - 用于局部波动性估计
    Rolling_MAD_D = calculate_sparse_rolling_mad(Z_Composite_D, N_WINDOW, K_STEP),
    Rolling_MAD_O18 = calculate_sparse_rolling_mad(Z_Composite_O18, N_WINDOW, K_STEP),
    
    # 动态控制限 (CL)
    CL_D = Rolling_Median_D + K_CL_SIGMA * Rolling_MAD_D,
    CL_O18 = Rolling_Median_O18 + K_CL_SIGMA * Rolling_MAD_O18
  )

# --- 3.2 Outlier 判定逻辑 ---
# 判定逻辑：Outlier = 超过 SL (静态规格限) OR 超过 CL (动态控制限)

df_qc_a <- df_qc_a %>%
  mutate(
    # (1) 静态规格限 (SL) 判定 (Z-Score尺度上的 SL 就是 K_SL_SIGMA)
    SL_D_Outlier = Z_Composite_D > K_SL_SIGMA,
    SL_O18_Outlier = Z_Composite_O18 > K_SL_SIGMA,
    
    # (2) 动态控制限 (CL) 判定 (仅对 Cycle_Seq >= N_WINDOW 的数据进行判定)
    CL_D_Outlier = (Z_Composite_D > CL_D) & (Cycle_Seq >= N_WINDOW),
    CL_O18_Outlier = (Z_Composite_O18 > CL_O18) & (Cycle_Seq >= N_WINDOW),
    
    # (3) 最终 QC-A Outlier 判定 (Total_Outlier = SL OR CL)
    Total_D_Outlier = SL_D_Outlier | CL_D_Outlier,
    Total_O18_Outlier = SL_O18_Outlier | CL_O18_Outlier,
    
    # Cycle 级别的最终判定：D 或 O18 任一异常，则该 Cycle 异常
    QC_A_Cycle_Outlier = Total_D_Outlier | Total_O18_Outlier
  )


# ======================================================================
# 4. 输出结果与统计摘要
# ======================================================================

str(df_qc_a)
table(df_qc_a$CL_D_Outlier)
message("\n--- 阶段 I：QC-A 质控结果 (部分展示) ---")
print(head(df_qc_a %>% 
             select(group.id, year, Cycle_Seq, Z_Composite_D, CL_D, K_SL_SIGMA, Total_D_Outlier, QC_A_Cycle_Outlier), 
           20))

message(paste0("\n--- 质控参数摘要：N=", N_WINDOW, ", K_CL=", K_CL_SIGMA, ", K_STEP=", K_STEP, ", K_SL=", K_SL_SIGMA, " ---"))

# 统计摘要
df_summary <- df_qc_a %>%
  summarise(
    Total_Cycles = n(),
    D_Outlier_Count = sum(Total_D_Outlier, na.rm = TRUE),
    O18_Outlier_Count = sum(Total_O18_Outlier, na.rm = TRUE),
    Combined_Outlier_Count = sum(QC_A_Cycle_Outlier, na.rm = TRUE),
    
    D_Outlier_Rate = round(D_Outlier_Count / Total_Cycles * 100, 2),
    O18_Outlier_Rate = round(O18_Outlier_Count / Total_Cycles * 100, 2),
    Combined_Outlier_Rate = round(Combined_Outlier_Count / Total_Cycles * 100, 2)
  )

message("\n--- QC-A 异常值总览 (总 Cycles) ---")
print(df_summary)

message("\n--- 窗口期 CL 异常值统计 (仅在 Cycle_Seq >= N_WINDOW 的数据上进行 CL 判定) ---")
df_cl_summary <- df_qc_a %>%
  filter(Cycle_Seq >= N_WINDOW) %>%
  summarise(
    Total_CL_Testable_Cycles = n(),
    CL_D_Count = sum(CL_D_Outlier, na.rm = TRUE),
    CL_O18_Count = sum(CL_O18_Outlier, na.rm = TRUE),
    CL_Combined_Count = sum(CL_D_Outlier | CL_O18_Outlier, na.rm = TRUE),
    
    CL_Combined_Rate = round(CL_Combined_Count / Total_CL_Testable_Cycles * 100, 2)
  )

print(df_cl_summary)



####### plot 

# --- 0. 确保加载所需的库 (如果尚未加载) ---
if (!requireNamespace("latex2exp", quietly = TRUE)) install.packages("latex2exp")
library(latex2exp)

# --- 1. 定义增强的 Base R 绘图函数 (英文标签) ---
# 绘制 Z_Composite 的时间序列，并标注 CL 和 SL 阈值
plot_z_time_series_base_r <- function(data, Z_col, CL_col, SL_col,
                                      K_SL_SIGMA, K_CL_SIGMA, Y_label_expr, Title_Suffix) {
  
  # Data extraction
  Z_data <- data[[Z_col]]
  CL_data <- data[[CL_col]]
  SL_threshold_value <- unique(data[[SL_col]])[1]
  Cycle_Seq <- data$Cycle_Seq
  
  # Outlier detection logic: Exceeding CL OR Exceeding SL
  is_outlier <- Z_data > CL_data | Z_data > SL_threshold_value
  is_normal <- !is_outlier
  
  # Set up plot margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # 1. Initialize the plot (Z_Composite time series)
  plot(Cycle_Seq, Z_data, 
       type = "n", # Do not plot yet, set up the canvas
       xlab = "Cycle Sequence Number (Time)",
       ylab = "", # Set later with title() for math expression
       main = paste0("Cycle-Level Z-Composite Time Series: ", Title_Suffix),
       ylim = range(c(Z_data, CL_data, SL_threshold_value), na.rm = TRUE),
       las = 1, # Horizontal axis labels
       cex.axis = 0.8,
       cex.lab = 1.2,
       cex.main = 1.2)
  
  # Add grid lines for readability
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 1, equilogs = TRUE)
  
  # Add Y-axis label with expression
  title(ylab = Y_label_expr, line = 3.5, cex.lab = 1.2)
  
  # Add the main connecting line
  lines(Cycle_Seq, Z_data, col = "gray50", lwd = 0.5)
  
  # 2. Plot Normal Cycles (Non-Outliers) - NEW COLORING
  points(Cycle_Seq[is_normal], Z_data[is_normal], 
         pch = 19, # Solid circle
         cex = 0.5, 
         col = "darkblue") # Normal points in dark blue
  
  # 3. Plot Outlier Cycles (Exceeding CL or SL) - ENHANCED HIGHLIGHT
  points(Cycle_Seq[is_outlier], Z_data[is_outlier], 
         pch = 21, # Circle with fill
         bg = "red", # Fill color: Outliers in RED
         col = "black", # Border color
         cex = 1.2) # Clearly visible size
  
  # 4. Add Dynamic Control Limit (CL) - Blue Dashed
  lines(Cycle_Seq, CL_data, 
        col = "blue", 
        lty = 2, # Dashed line
        lwd = 1.5)
  
  # 5. Add Static Specification Limit (SL) - Red Solid
  abline(h = SL_threshold_value, 
         col = "red", 
         lty = 1, # Solid line
         lwd = 2)
  
  # 6. Annotations (Control Limit Equations)
  max_seq <- max(Cycle_Seq, na.rm = TRUE)
  
  # CL Annotation
  text(x = max_seq * 0.95, 
       y = CL_data[length(CL_data[!is.na(CL_data)])] + 0.02 * diff(par("usr")[3:4]), # Offset from the end of the line
       labels = bquote(CL: ~ Median["Window"] + .(K_CL_SIGMA) * MAD["Window"]),
       col = "blue", 
       adj = c(1, 0), # Right-align, bottom-align
       cex = 0.8)
  
  # SL Annotation
  text(x = max_seq * 0.95, 
       y = SL_threshold_value + 0.06 * diff(par("usr")[3:4]), # Offset slightly from the line
       labels = bquote(SL: ~ Median["Global"] + .(K_SL_SIGMA) * MAD["Global"]),
       col = "red", 
       adj = c(1, 0), # Right-align, bottom-align
       cex = 0.8)
  
  # Title Subtitle (Static/Dynamic Thresholds)
  mtext(paste0("CL (Dynamic): Rolling Median + ", K_CL_SIGMA, "*MAD | SL (Static): Global Median + ", K_SL_SIGMA, "*MAD"), 
        side = 3, line = 0.5, cex = 0.9)

  # 7. Add Legend - NEW
  legend_text <- c("Normal Cycle", "Outlier Cycle", "Dynamic CL", "Static SL")
  legend_colors <- c("darkblue", "black", "blue", "red")
  legend_pch <- c(19, 21, NA, NA)
  legend_lty <- c(NA, NA, 2, 1)
  legend_lwd <- c(NA, NA, 1.5, 2)
  legend_pt_bg <- c(NA, "red", NA, NA) # Fill color for the outlier circle
  
  legend("topleft", 
         legend = legend_text, 
         col = legend_colors, 
         lty = legend_lty, 
         lwd = legend_lwd,
         pch = legend_pch,
         pt.bg = legend_pt_bg,
         bty = "n", # Remove box around legend
         cex = 0.9)
  
  # Reset margins to default (optional)
  par(mar = c(5, 4, 4, 2) + 0.1)
}


# --- 2. 假设 df_qc_a 数据框已存在 (您的 QC-A 脚本执行结果) ---
# **重要：请确保在运行此绘图脚本之前，您已成功运行 QC-A 脚本并生成了 df_qc_a**
# 为了让这段代码能够运行，我保留了原有的结构。
K_CL_SIGMA_VALUE <- 3 # 保持您的默认值

if (exists("df_qc_a")) {
  
  # 提取 K_SL_SIGMA (如果存在)
  K_SL_SIGMA_VALUE <- unique(df_qc_a$K_SL_SIGMA)[1]
  if (is.null(K_SL_SIGMA_VALUE) || is.na(K_SL_SIGMA_VALUE)) {
    # Fallback to default if not found in data frame
    K_SL_SIGMA_VALUE <- 5 
  }
  
  OUTPUT_DIR <- "./WVISS_QC_Plots"
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
  
  # 3.1 绘制 dD 的 Z-Composite 时间序列
  pdf(file.path(OUTPUT_DIR, "Z_Composite_dD_TimeSeries_BaseR_Enhanced.pdf"), width = 14, height = 7)
  plot_z_time_series_base_r(
    data = df_qc_a,
    Z_col = "Z_Composite_D",
    CL_col = "CL_D", 
    SL_col = "SL_D", 
    K_SL_SIGMA = K_SL_SIGMA_VALUE,
    K_CL_SIGMA = K_CL_SIGMA_VALUE,
    # Use expression for the Y-axis label (dD)
    Y_label_expr = expression(Z["Composite"] * " (" * delta * D * ")"),
    Title_Suffix = "Hydrogen Isotope (dD)"
  )
  dev.off() # Close the device to save the file
  message("Z_Composite_dD_TimeSeries_BaseR_Enhanced.pdf Plotting completed.")
  
  # 3.2 绘制 d18O 的 Z-Composite 时间序列
  pdf(file.path(OUTPUT_DIR, "Z_Composite_d18O_TimeSeries_BaseR_Enhanced.pdf"), width = 14, height = 7)
  plot_z_time_series_base_r(
    data = df_qc_a,
    Z_col = "Z_Composite_O18",
    CL_col = "CL_O18", 
    SL_col = "SL_O18", 
    K_SL_SIGMA = K_SL_SIGMA_VALUE,
    K_CL_SIGMA = K_CL_SIGMA_VALUE,
    # Use expression for the Y-axis label (d18O)
    Y_label_expr = expression(Z["Composite"] * " (" * delta^18 * O * ")"),
    Title_Suffix = "Oxygen Isotope (d18O)"
  )
  dev.off() # Close the device to save the file
  message("Z_Composite_d18O_TimeSeries_BaseR_Enhanced.pdf Plotting completed.")
  
} else {
  message("Error: Data frame 'df_qc_a' not found. Please run the QC-A script first to generate the QC results data.")
}
