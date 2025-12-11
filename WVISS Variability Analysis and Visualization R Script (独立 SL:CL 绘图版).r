# ======================================================================
# WVISS Variability Analysis and Visualization R Script (独立 SL/CL 绘图版)
# ======================================================================

# ----------------------------------------------------------------------
# 0. 环境设置与库加载
# ----------------------------------------------------------------------
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
if (!requireNamespace("maditr", quietly = TRUE)) install.packages("maditr")
if (!requireNamespace("latex2exp", quietly = TRUE)) install.packages("latex2exp")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(tidyverse)
library(zoo)
library(maditr) # 用于 dt_mutate
library(latex2exp)
library(scales) 
tryCatch(library(RcppRoll), error = function(e) message("RcppRoll not found. Using zoo::rollapply (slower) for rolling median."))

options("digits" = 9) 


# 设置输出路径
OUTPUT_DIR <- "./WVISS_QC_Results"
VISUALIZATION_DIR <- file.path(OUTPUT_DIR, "Visualizations")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
if (!dir.exists(VISUALIZATION_DIR)) dir.create(VISUALIZATION_DIR)

# ----------------------------------------------------------------------
# 1. 数据加载与预处理 (沿用用户提供的代码结构)
# ----------------------------------------------------------------------

# ACTION REQUIRED: Manually set the path where WVISS.all.csv is located
wd.synthesis <- '/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'
data_path <- file.path(wd.synthesis, 'WVISS.all.csv')

wd.QC_A <- '/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/'

# 更改工作目录到数据所在路径
setwd(wd.synthesis)

WVISS.all <- NULL
if (file.exists(data_path)) {
  WVISS.all <- read.csv(data_path, header = T, stringsAsFactors = FALSE)
  WVISS.all <- as_tibble(WVISS.all) 
  message("WVISS.all.csv data loaded successfully.")
} else {
  stop(paste("Error: WVISS.all.csv file not found at:", data_path, "\n"))
}

# --- Cycle Chronological Sorting Logic ---
WVISS.data <- WVISS.all %>%
  # 1. Calculate the mean time (DOY) for each group.id (Cycle)
  group_by(group.id) %>%
  mutate(group_DOY_mean = mean(DOY_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  
  # 2. Sort data by mean DOY to ensure correct chronological order
  arrange(group_DOY_mean)

# 3. Create Cycle_Seq (the chronological cycle number) and Date columns
chronological_groups <- unique(WVISS.data$group.id)

WVISS.data <- WVISS.data %>%
  mutate(
    Start_Date = as.Date("2011-01-01"),
    # Calculate Date based on DOY_mean
    Date = Start_Date + floor(DOY_mean) - 1,
    # Assign Cycle_Seq based on chronological order
    Cycle_Seq = as.numeric(factor(group.id, levels = chronological_groups)),
    # Extract Year and Month for aggregation (Crucial for Step 5)
    Month = as.numeric(format(Date, "%m")),
    Year = as.numeric(format(Date, "%Y"))
  ) %>%
  select(-Start_Date, -group_DOY_mean)

message("Data pre-processing and chronological sorting completed.")
setwd(wd.QC_A) # 更改工作目录到 QC scripts 所在路径
str(WVISS.data )

# ----------------------------------------------------------------------
# 2. QC 参数定义 (全局阈值 k, 滑动窗口参数 N, 步长 s)
# ----------------------------------------------------------------------

# 全局（SL - Static/System Level）和 局部（WL/CL - Window Level/Control Limit）离群值阈值
K_SL_SIGMA <- 5.0 # 全局离群值阈值的 Sigma 倍数 (SL 阈值: 5 倍)
K_WL_SIGMA <- 3.0 # 局部离群值阈值的 Sigma 倍数 (CL 阈值: 3 倍)

# 滑动窗口参数 (CL)
WINDOW_STEP <- 5L    # 滑动窗口步长 (Cycle 数)
WINDOW_WIDTH <- 270L # 滑动窗口宽度 (Cycle 数)

ISOTOPES <- c("D", "O18") # 同位素种类

message(paste("QC 参数设置：K_SL_SIGMA =", K_SL_SIGMA, 
              ", K_WL_SIGMA =", K_WL_SIGMA, 
              ", N (Window Width) =", WINDOW_WIDTH, 
              ", s (Step) =", WINDOW_STEP))


# ----------------------------------------------------------------------
# 3. Cycle 级别变异性（SD 和 Range）计算
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# 3. Cycle 级别变异性（SD 和 Range）计算 (【修正版】: 移除严格过滤)
# ----------------------------------------------------------------------
df_cycle_metrics <- WVISS.data %>%
  # 保持多列分组：确保 Date, Year, Month, Cycle_Seq 都被保留（元数据透传）
  group_by(group.id, Date, Year, Month, Cycle_Seq) %>% 
  summarise(
    # dD 变异性
    Ddel_mean = mean(Ddel_mean, na.rm = TRUE), # 保留 Cycle 均值
    D_SD = sd(Ddel_mean, na.rm = TRUE),
    D_Range = max(Ddel_mean, na.rm = TRUE) - min(Ddel_mean, na.rm = TRUE),
    # dO18 变异性
    O18del_mean = mean(O18del_mean, na.rm = TRUE), # 保留 Cycle 均值
    O18_SD = sd(O18del_mean, na.rm = TRUE),
    O18_Range = max(O18del_mean, na.rm = TRUE) - min(O18del_mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # 显式地按 Cycle_Seq 排序
  arrange(Cycle_Seq) 
  
  # *** 关键修正: 移除了 filter_at() 语句 ***
  # 之前移除所有 SD 或 Range 为 NA 的 Cycle，导致数据丢失。

message(paste("Cycle 级别 SD/Range 计算完成。【修正后】行数:", nrow(df_cycle_metrics), "个 Cycle。"))


# ======================================================================
# 4. 离群值识别与 Z-score 计算 【已修正 CL 缺失问题】
# ======================================================================

# --- A. 全局（SL）统计量：Median & MAD ---
# 计算 SD 和 Range 的全局中值（Median）和中值绝对偏差（MAD）
df_global_stats <- df_cycle_metrics %>%
  summarise(
    D_SD_Median = median(D_SD, na.rm = TRUE), D_SD_MAD = mad(D_SD, na.rm = TRUE),
    D_Range_Median = median(D_Range, na.rm = TRUE), D_Range_MAD = mad(D_Range, na.rm = TRUE),
    
    O18_SD_Median = median(O18_SD, na.rm = TRUE), O18_SD_MAD = mad(O18_SD, na.rm = TRUE),
    O18_Range_Median = median(O18_Range, na.rm = TRUE), O18_Range_MAD = mad(O18_Range, na.rm = TRUE)
  )

# --- B. 局部（WL/CL）鲁棒滑动 MAD 函数 (替换了原来的复杂 B 部分) ---
# 此函数用于 zoo::rollapply，确保 MAD 即使在 NA 很多或数据不活跃时也能鲁棒。
roll_mad <- function(x) {
    # 仅使用非 NA 值
    x_valid <- x[!is.na(x)]
    # 如果有效数据点太少 (例如少于 5 个)，则返回 NA
    if (length(x_valid) < 5) { 
        return(NA_real_)
    }
    # 计算 MAD (constant=1.4826 确保其与 SD 一致)
    mad(x_valid, center = median(x_valid), constant = 1.4826)
}


# --- C. 离群值识别和 Z-score 计算 (SL & CL 完整鲁棒版) ---
# 扩展数据框以包含全局统计量和 Z-score

df_metrics_extended <- df_cycle_metrics %>%
  # ----------------------------------------------------------------------
  # SL (全局) Z-Score 计算 - 【鲁棒修正版】
  # ----------------------------------------------------------------------
  mutate(
    # --- D 同位素 SL Z-Scores ---
    Z_D_Range_SL = case_when(
      # MAD > 0
      df_global_stats$D_Range_MAD[1] > 0 ~ (D_Range - df_global_stats$D_Range_Median[1]) / df_global_stats$D_Range_MAD[1],
      # MAD = 0 (所有非中位数的值视为 Inf 离群)
      df_global_stats$D_Range_MAD[1] == 0 ~ ifelse(D_Range == df_global_stats$D_Range_Median[1], 0, Inf),
      TRUE ~ NA_real_
    ),
    # D_SD SL Z-score (可能 MAD=NA)
    Z_D_SD_SL = (D_SD - df_global_stats$D_SD_Median[1]) / df_global_stats$D_SD_MAD[1],
    
    # --- O18 同位素 SL Z-Scores ---
    Z_O18_Range_SL = case_when(
      df_global_stats$O18_Range_MAD[1] > 0 ~ (O18_Range - df_global_stats$O18_Range_Median[1]) / df_global_stats$O18_Range_MAD[1],
      df_global_stats$O18_Range_MAD[1] == 0 ~ ifelse(O18_Range == df_global_stats$O18_Range_Median[1], 0, Inf),
      TRUE ~ NA_real_
    ),
    Z_O18_SD_SL = (O18_SD - df_global_stats$O18_SD_Median[1]) / df_global_stats$O18_SD_MAD[1],
    
    # SL Composite (取 Range 和 SD 中绝对值的最大值)
    Z_Composite_D_SL = pmax(abs(Z_D_Range_SL), abs(Z_D_SD_SL), na.rm = TRUE),
    Z_Composite_O18_SL = pmax(abs(Z_O18_Range_SL), abs(Z_O18_SD_SL), na.rm = TRUE),
    
    # 替换 Inf 为 K_SL_SIGMA + 1
    Z_Composite_D_SL = ifelse(is.infinite(Z_Composite_D_SL), K_SL_SIGMA + 1, Z_Composite_D_SL),
    Z_Composite_O18_SL = ifelse(is.infinite(Z_Composite_O18_SL), K_SL_SIGMA + 1, Z_Composite_O18_SL),
    
    # SL 离群值标记
    Outlier_D_SL = as.numeric(!is.na(Z_Composite_D_SL) & Z_Composite_D_SL > K_SL_SIGMA),
    Outlier_O18_SL = as.numeric(!is.na(Z_Composite_O18_SL) & Z_Composite_O18_SL > K_SL_SIGMA)
  ) %>%
  
  # ----------------------------------------------------------------------
  # CL (滑动窗口) 统计量和 Z-Score 计算 - 【缺失补充】
  # ----------------------------------------------------------------------
  # 1. 计算滑动窗口统计量 (使用 trailing window: align = "right")
  # 注意：SD 滑动统计量因 Cycle 变异性数据结构问题可能大部分为 NA，故此处仅计算 Range
  mutate(
    # D Range
    D_Range_Median_CL = zoo::rollapply(D_Range, width = WINDOW_WIDTH, FUN = median, na.rm = TRUE, align = "right", fill = NA),
    D_Range_MAD_CL = zoo::rollapply(D_Range, width = WINDOW_WIDTH, FUN = roll_mad, align = "right", fill = NA),
    # O18 Range
    O18_Range_Median_CL = zoo::rollapply(O18_Range, width = WINDOW_WIDTH, FUN = median, na.rm = TRUE, align = "right", fill = NA),
    O18_Range_MAD_CL = zoo::rollapply(O18_Range, width = WINDOW_WIDTH, FUN = roll_mad, align = "right", fill = NA)
  ) %>%
  
  # 2. 计算 CL 模式 Z-score (处理 MAD=0 的鲁棒逻辑)
  mutate(
    # --- D 同位素 CL Z-Scores (基于 Range) ---
    Z_D_Range_CL = case_when(
      D_Range_MAD_CL > 0 ~ (D_Range - D_Range_Median_CL) / D_Range_MAD_CL,
      D_Range_MAD_CL == 0 ~ ifelse(D_Range == D_Range_Median_CL, 0, Inf),
      TRUE ~ NA_real_
    ),
    
    # --- O18 同位素 CL Z-Scores (基于 Range) ---
    Z_O18_Range_CL = case_when(
      O18_Range_MAD_CL > 0 ~ (O18_Range - O18_Range_Median_CL) / O18_Range_MAD_CL,
      O18_Range_MAD_CL == 0 ~ ifelse(O18_Range == O18_Range_Median_CL, 0, Inf),
      TRUE ~ NA_real_
    ),
    
    # 3. 计算 CL Composite Z-score (仅使用 Range 组件)
    Z_Composite_D_CL = abs(Z_D_Range_CL), 
    Z_Composite_O18_CL = abs(Z_O18_Range_CL),
    
    # 4. 替换 Inf 值 (为了 Outlier 判断和绘图稳定)
    Z_Composite_D_CL = ifelse(is.infinite(Z_Composite_D_CL), K_WL_SIGMA + 1, Z_Composite_D_CL),
    Z_Composite_O18_CL = ifelse(is.infinite(Z_Composite_O18_CL), K_WL_SIGMA + 1, Z_Composite_O18_CL),
    
    # 5. CL 离群值标记
    Outlier_D_CL = as.numeric(!is.na(Z_Composite_D_CL) & Z_Composite_D_CL > K_WL_SIGMA),
    Outlier_O18_CL = as.numeric(!is.na(Z_Composite_O18_CL) & Z_Composite_O18_CL > K_WL_SIGMA)
  )

# ======================================================================
# 5. 诊断 Step 4 结果 (确保 CL 列已创建)
# ======================================================================
message("\n--- Step 4 结果诊断 (验证 CL 列) ---")
if (exists("df_metrics_extended") && nrow(df_metrics_extended) > 0) {
    total_cycles <- nrow(df_metrics_extended)
    na_sl_d <- sum(is.na(df_metrics_extended$Z_Composite_D_SL))
    na_cl_d <- sum(is.na(df_metrics_extended$Z_Composite_D_CL))
    
    message(paste("总 Cycle 数:", total_cycles))
    message(paste("Z_Composite_D_SL 列 NA 比例:", round(na_sl_d / total_cycles * 100, 2), "%"))
    message(paste("Z_Composite_D_CL 列 NA 比例 (预期前 N-1 个 NA):", round(na_cl_d / total_cycles * 100, 2), "%"))

    if (na_sl_d == total_cycles) {
        message("SL 结论：所有 SL Z-score 仍为 NA。请检查 df_global_stats 的 D_SD/O18_SD 是否为 NA。")
    }
    if (na_cl_d == total_cycles) {
        message("CL 结论：所有 CL Z-score 均为 NA。请检查 WINDOW_WIDTH 是否太大或 D_Range 是否全为 NA。")
    }
} else {
    message("警告：df_metrics_extended 未找到或为空。")
}


# ----------------------------------------------------------------------
# 6. 时间序列图绘制（SL 和 CL 独立趋势）- 【Base R 基础绘图版】
# ----------------------------------------------------------------------

# 绘图函数：用于统一绘制 SL 或 CL Z-score 时间序列图 (使用 Base R)
# mode: "SL" 或 "CL"
plot_z_time_series_base_r <- function(data, z_col, plot_k_sigma, y_label, title_label, mode) {
  
  # 1. 检查数据
  if (nrow(data) == 0 || all(is.na(data$Date))) {
      warning(paste("Plot Failed for", title_label, mode, ": Input data frame is empty or contains no valid 'Date' values."))
      return(invisible(NULL)) 
  }
  
  # 2. 准备数据和参数
  Z_threshold <- plot_k_sigma 
  Outlier_Status_Col <- paste0("Outlier_", ifelse(grepl("_D_", z_col), "D", "O18"), "_", mode)
  
  # 过滤掉 NA Z-score 的行（Base R 绘图的关键步骤）
  plot_data <- data[!is.na(data[[z_col]]), ]
  
  # 检查 Outlier Status 列是否存在
  if (!Outlier_Status_Col %in% names(data)) {
      warning(paste("Plot Failed for", title_label, mode, ": Missing outlier status column:", Outlier_Status_Col))
      return(invisible(NULL))
  }
  
  if (nrow(plot_data) == 0) {
     warning(paste("Plot Failed for", title_label, mode, ": No non-NA Z-scores found after filtering."))
     return(invisible(NULL))
  }
  
  # 颜色映射: 0=Blue (Normal), 1=Red (Outlier)
  outlier_colors <- c("0" = "blue", "1" = "red")
  plot_colors <- outlier_colors[as.character(plot_data[[Outlier_Status_Col]])]
  
  # 3. 初始化绘图窗口
  
  # 标题/副标题设置
  main_title <- paste(title_label, mode, "时间序列")
  sub_text <- paste0(mode, "模式：", ifelse(mode == "CL", 
                                            paste0("滑动窗口 N=", WINDOW_WIDTH, ", s=", WINDOW_STEP, "；"), 
                                            "全局静态统计；"), 
                     "阈值 ", plot_k_sigma, "\u03c3") # \u03c3 是 unicode sigma
  
  # 绘制所有 Z-score 的范围
  y_range <- range(c(0, Z_threshold, plot_data[[z_col]]), na.rm = TRUE)
  
  plot(plot_data$Date, plot_data[[z_col]], 
       type = "n", # 不绘制任何点，只设置坐标轴
       xlab = "日期", 
       ylab = y_label, 
       main = main_title,
       sub = sub_text,
       ylim = y_range,
       xaxt = "n" # 抑制默认 X 轴，使用 Date 格式
  )
  
  # 4. 绘制时间序列
  
  # 绘制连接线 (使用中性蓝色)
  lines(plot_data$Date, plot_data[[z_col]], col = "blue", lty = 1, lwd = 0.5) 
  
  # 绘制点，颜色根据离群值状态决定
  points(plot_data$Date, plot_data[[z_col]], 
         col = plot_colors, 
         pch = 16, # 实心圆
         cex = 0.8  # 点的大小
  )
  
  # 5. X 轴日期格式化
  axis.Date(1, plot_data$Date, format = "%Y-%m", cex.axis = 0.8)

  # 6. 阈值横线 (`abline` 替换 `geom_hline`)
  abline(h = Z_threshold, lty = 2, col = "red", lwd = 1)
  
  # 7. 阈值文字 (`text` 替换 `geom_text`)
  Date_max <- max(plot_data$Date, na.rm = TRUE)
  # 手动计算 X 轴位置，避免使用 text_offset 这种不可靠的方法
  x_pos_text <- max(plot_data$Date, na.rm = TRUE)
  y_pos_text <- Z_threshold
  
  threshold_label <- paste0(mode, " Threshold: ", plot_k_sigma, "\u03c3")
  
  text(x = x_pos_text, y = y_pos_text, 
       labels = threshold_label, 
       adj = c(1, -0.5), # 右对齐, 略微上浮
       col = "red", 
       cex = 0.8
  )
  
  # 8. 图例 (`legend` 替换 `scale_color_manual` 的图例)
  legend("bottomright", 
         legend = c("Normal Cycle", "Outlier Cycle"), 
         col = c("blue", "red"), 
         pch = 16, 
         lty = 1, 
         title = paste0(mode, " Status"), 
         bty = "n", # 无边框
         cex = 0.9
  )
}

# ----------------------------------------------------------------------
# 6. (续) 调用绘图函数和保存文件
# ----------------------------------------------------------------------

# --- 1. 绘制 SL (全局) 时间序列 (阈值 K_SL_SIGMA = 5.0) ---
str(df_metrics_extended)
# Z_Composite_D (SL)
pdf(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_D_SL_Separate_BaseR.pdf"), width = 12, height = 6)
plot_z_time_series_base_r(df_metrics_extended, "Z_Composite_D_SL", K_SL_SIGMA, 
                          "Z-Score (dD)", "Z-Score of dD Composite Variability", "SL")
dev.off() # 关闭图形设备，保存文件
message("WVISS_Z_Composite_D_SL_Separate_BaseR.pdf 绘制完成。")

# Z_Composite_O18 (SL)
pdf(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_O18_SL_Separate_BaseR.pdf"), width = 12, height = 6)
plot_z_time_series_base_r(df_metrics_extended, "Z_Composite_O18_SL", K_SL_SIGMA, 
                          "Z-Score (d18O)", "Z-Score of d18O Composite Variability", "SL")
dev.off() # 关闭图形设备，保存文件
message("WVISS_Z_Composite_O18_SL_Separate_BaseR.pdf 绘制完成。")


# --- 2. 绘制 CL (滑动窗口) 时间序列 (阈值 K_WL_SIGMA = 3.0) ---

# Z_Composite_D (CL)
pdf(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_D_CL_Separate_BaseR.pdf"), width = 12, height = 6)
plot_z_time_series_base_r(df_metrics_extended, "Z_Composite_D_CL", K_WL_SIGMA, 
                          "Z-Score (dD)", "Z-Score of dD Composite Variability", "CL")
dev.off() # 关闭图形设备，保存文件
message("WVISS_Z_Composite_D_CL_Separate_BaseR.pdf 绘制完成。")

# Z_Composite_O18 (CL)
pdf(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_O18_CL_Separate_BaseR.pdf"), width = 12, height = 6)
plot_z_time_series_base_r(df_metrics_extended, "Z_Composite_O18_CL", K_WL_SIGMA, 
                          "Z-Score (d18O)", "Z-Score of d18O Composite Variability", "CL")
dev.off() # 关闭图形设备，保存文件
message("WVISS_Z_Composite_O18_CL_Separate_BaseR.pdf 绘制完成。")


message("\n所有 4 个 SL/CL 独立 Z-score 时间序列图已生成 (Base R 版本)。")