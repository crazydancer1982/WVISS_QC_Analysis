# WVISS Variability Analysis and Visualization R Script (Final Optimized Version).r
#
# 作者：[您的姓名/组织]
# 日期：2025-12-09
# 目标：计算 WVISS 数据的周期（Cycle）级别变异性（SD/Range），
#      识别滑动窗口和全局离群值，计算综合 Z-score，并进行可视化。
# 更新内容：
# 1. 调整参数：K_SL_SIGMA = 5.0，WINDOW_STEP = 5，WINDOW_WIDTH = 270。
# 2. 时间序列图：增加 Z_Composite_O18，阈值横线添加 Median+5*MAD 注释。
# 3. 堆叠图：增加 Z_Composite_O18 和 Z_Composite_Combined (D或O18超限) 两个堆叠图。
#
# ----------------------------------------------------------------------
# 0. 环境设置与库加载
# ----------------------------------------------------------------------
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
if (!requireNamespace("maditr", quietly = TRUE)) install.packages("maditr")
if (!requireNamespace("latex2exp", quietly = TRUE)) install.packages("latex2exp")

library(tidyverse)
library(zoo)
library(maditr)
library(latex2exp)

# 设置工作目录和数据路径
# 请根据您的实际情况修改 DATA_PATH
DATA_PATH <- "WVISS_Raw_Data.csv" # 假设数据文件名为 WVISS_Raw_Data.csv
OUTPUT_DIR <- "WVISS_QC_Results"
VISUALIZATION_DIR <- file.path(OUTPUT_DIR, "Visualizations")

# 确保输出目录存在
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
if (!dir.exists(VISUALIZATION_DIR)) dir.create(VISUALIZATION_DIR)

# ----------------------------------------------------------------------
# 1. QC 参数定义 (全局阈值 k, 滑动窗口参数 N, 步长 s)
# ----------------------------------------------------------------------
# 全局（SL - Static/System Level）和 局部（WL - Window Level）离群值阈值
K_SL_SIGMA <- 5.0 # 全局离群值阈值的 Sigma 倍数 (更新为 5 倍)
K_WL_SIGMA <- 3.0 # 局部离群值阈值的 Sigma 倍数 (保持 3 倍)

# 滑动窗口参数 (WL)
WINDOW_STEP <- 5L    # 滑动窗口步长 (Cycle 数) (更新为 5 Cycles)
WINDOW_WIDTH <- 270L # 滑动窗口宽度 (Cycle 数) (更新为 270 Cycles)

message(paste("QC 参数设置：K_SL_SIGMA =", K_SL_SIGMA, 
              ", K_WL_SIGMA =", K_WL_SIGMA, 
              ", N (Window Width) =", WINDOW_WIDTH, 
              ", s (Step) =", WINDOW_STEP))

# ----------------------------------------------------------------------
# 2. 数据读取与预处理
# ----------------------------------------------------------------------
if (!file.exists(DATA_PATH)) {
  stop(paste("未找到数据文件：", DATA_PATH, "请检查路径是否正确。"))
}

# 假设数据包含以下关键列：
# 'Date' (日期时间), 'Year', 'Month', 'group.id' (Cycle ID), 
# 'Ddel_mean' (dD 均值), 'O18del_mean' (dO18 均值)
WVISS.data <- read.csv(DATA_PATH) %>%
  mutate(
    Date = as.Date(Date), # 确保 Date 为 Date 类型
    # 假设 'group.id' 已经定义好
    group.id = as.character(group.id) 
  )

# 按时间顺序（假设 Date/Year/Month/DOY）生成 Cycle 序列号
# WVISS.data 应该在之前的步骤中已经加入了 DOY 等列用于排序
# 为确保鲁棒性，这里假设 WVISS.data 已经按时间排序
WVISS.data <- WVISS.data %>%
  arrange(Date, Year, Month, group.id) %>%
  # 基于 group.id 创建唯一的 Cycle 序列号
  mutate(Cycle_Seq = as.integer(factor(group.id, levels = unique(group.id)))) 

message(paste("数据读取完成。共计", nrow(WVISS.data), "行数据。"))

# ----------------------------------------------------------------------
# 3. Cycle 级别变异性（SD 和 Range）计算
# ----------------------------------------------------------------------
df_cycle_metrics <- WVISS.data %>%
  # 保持多列分组：确保 Date, Year, Month, Cycle_Seq 都被保留（元数据透传）
  group_by(group.id, Date, Year, Month, Cycle_Seq) %>% 
  summarise(
    # dD 变异性
    D_SD = sd(Ddel_mean, na.rm = TRUE),
    D_Range = max(Ddel_mean, na.rm = TRUE) - min(Ddel_mean, na.rm = TRUE),
    # dO18 变异性
    O18_SD = sd(O18del_mean, na.rm = TRUE),
    O18_Range = max(O18del_mean, na.rm = TRUE) - min(O18del_mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # 显式地按 Cycle_Seq 排序
  arrange(Cycle_Seq) %>%
  # 移除 SD 或 Range 为 Inf 的行（通常发生在 Cycle 只有一个数据点时）
  filter_at(vars(D_SD, D_Range, O18_SD, O18_Range), all_vars(!is.infinite(.)))

message(paste("Cycle 级别 SD/Range 计算完成。共计", nrow(df_cycle_metrics), "个 Cycle。"))

# ----------------------------------------------------------------------
# 4. 离群值识别与 Z-score 计算
# ----------------------------------------------------------------------

# --- A. 全局（SL）统计量：Median & MAD ---
# 计算 SD 和 Range 的全局中值（Median）和中值绝对偏差（MAD）
df_global_stats <- df_cycle_metrics %>%
  summarise(
    D_SD_Median = median(D_SD, na.rm = TRUE),
    D_SD_MAD = mad(D_SD, na.rm = TRUE),
    D_Range_Median = median(D_Range, na.rm = TRUE),
    D_Range_MAD = mad(D_Range, na.rm = TRUE),
    
    O18_SD_Median = median(O18_SD, na.rm = TRUE),
    O18_SD_MAD = mad(O18_SD, na.rm = TRUE),
    O18_Range_Median = median(O18_Range, na.rm = TRUE),
    O18_Range_MAD = mad(O18_Range, na.rm = TRUE)
  )

# --- B. 局部（WL）滑动窗口统计量 ---
# 使用 zoo::rollapplyR 进行滑动窗口计算
# 注意：数据必须先按 Cycle_Seq 排序
roll_median_mad <- function(x) {
  # 计算 median 和 MAD，返回一个命名向量
  c(Median = median(x, na.rm = TRUE), MAD = mad(x, na.rm = TRUE))
}

# 计算滑动窗口统计量，并将结果添加到 df_cycle_metrics
df_cycle_metrics_wl <- df_cycle_metrics %>%
  dt_mutate(
    # dD_SD 局部统计量
    WL_D_SD_Median = rollapplyr(D_SD, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['Median'], by = WINDOW_STEP, fill = NA, align = "right"),
    WL_D_SD_MAD = rollapplyr(D_SD, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['MAD'], by = WINDOW_STEP, fill = NA, align = "right"),
    # dD_Range 局部统计量
    WL_D_Range_Median = rollapplyr(D_Range, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['Median'], by = WINDOW_STEP, fill = NA, align = "right"),
    WL_D_Range_MAD = rollapplyr(D_Range, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['MAD'], by = WINDOW_STEP, fill = NA, align = "right"),
    
    # dO18_SD 局部统计量
    WL_O18_SD_Median = rollapplyr(O18_SD, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['Median'], by = WINDOW_STEP, fill = NA, align = "right"),
    WL_O18_SD_MAD = rollapplyr(O18_SD, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['MAD'], by = WINDOW_STEP, fill = NA, align = "right"),
    # dO18_Range 局部统计量
    WL_O18_Range_Median = rollapplyr(O18_Range, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['Median'], by = WINDOW_STEP, fill = NA, align = "right"),
    WL_O18_Range_MAD = rollapplyr(O18_Range, width = WINDOW_WIDTH, FUN = function(x) roll_median_mad(x)['MAD'], by = WINDOW_STEP, fill = NA, align = "right")
  )

# C. 离群值识别和 Z-score 计算
# 扩展数据框以包含全局统计量
df_metrics_extended <- df_cycle_metrics_wl %>%
  dt_mutate(
    # 从全局统计量数据框中获取中值和 MAD
    D_SD_Median_SL = df_global_stats$D_SD_Median, D_SD_MAD_SL = df_global_stats$D_SD_MAD,
    D_Range_Median_SL = df_global_stats$D_Range_Median, D_Range_MAD_SL = df_global_stats$D_Range_MAD,
    O18_SD_Median_SL = df_global_stats$O18_SD_Median, O18_SD_MAD_SL = df_global_stats$O18_SD_MAD,
    O18_Range_Median_SL = df_global_stats$O18_Range_Median, O18_Range_MAD_SL = df_global_stats$O18_Range_MAD
  ) %>%
  dt_mutate(
    # 1. 计算 Z-score (SL)
    Z_D_SD_SL = (D_SD - D_SD_Median_SL) / D_SD_MAD_SL,
    Z_D_Range_SL = (D_Range - D_Range_Median_SL) / D_Range_MAD_SL,
    Z_O18_SD_SL = (O18_SD - O18_SD_Median_SL) / O18_SD_MAD_SL,
    Z_O18_Range_SL = (O18_Range - O18_Range_Median_SL) / O18_Range_MAD_SL,
    
    # 2. 计算 Z-score (WL)
    Z_D_SD_WL = (D_SD - WL_D_SD_Median) / WL_D_SD_MAD,
    Z_D_Range_WL = (D_Range - WL_D_Range_Median) / WL_D_Range_MAD,
    Z_O18_SD_WL = (O18_SD - WL_O18_SD_Median) / WL_O18_SD_MAD,
    Z_O18_Range_WL = (O18_Range - WL_O18_Range_Median) / WL_O18_Range_MAD
  ) %>%
  dt_mutate(
    # 3. 综合 Z-score (取 SL 和 WL Z-score 的较大值)
    Z_D_SD = pmax(Z_D_SD_SL, Z_D_SD_WL, na.rm = TRUE),
    Z_D_Range = pmax(Z_D_Range_SL, Z_D_Range_WL, na.rm = TRUE),
    Z_O18_SD = pmax(Z_O18_SD_SL, Z_O18_SD_WL, na.rm = TRUE),
    Z_O18_Range = pmax(Z_O18_Range_SL, Z_O18_Range_WL, na.rm = TRUE)
  ) %>%
  dt_mutate(
    # 4. 最终综合 Z-score (取 SD 和 Range 综合 Z-score 的较大值)
    # 这将用于最终的 QC 标记
    Z_Composite_D = pmax(Z_D_SD, Z_D_Range, na.rm = TRUE),
    Z_Composite_O18 = pmax(Z_O18_SD, Z_O18_Range, na.rm = TRUE)
  ) %>%
  dt_mutate(
    # 5. 确定全局阈值 (Median + k*MAD)
    # 注意：这里计算的是值域（Value Domain）的阈值
    D_Threshold_SL = D_SD_Median_SL + K_SL_SIGMA * D_SD_MAD_SL, # SD 的全局阈值，用于 dD
    O18_Threshold_SL = O18_SD_Median_SL + K_SL_SIGMA * O18_SD_MAD_SL, # SD 的全局阈值，用于 d18O
    
    # Z-score 的全局阈值，用于 Z-Composite 的判断
    Z_Threshold_SL = K_SL_SIGMA # Z-score 超过 K_SL_SIGMA 即为离群值
  ) %>%
  dt_mutate(
    # 6. 离群值标记
    # 使用 Z_Composite 超过 Z_Threshold_SL 来标记离群值 (1)
    Outlier_D = as.integer(Z_Composite_D > Z_Threshold_SL),
    Outlier_O18 = as.integer(Z_Composite_O18 > Z_Threshold_SL)
  )

message("离群值识别和 Z-score 计算完成。")

# ----------------------------------------------------------------------
# 5. 月度聚合与离群值百分比计算
# ----------------------------------------------------------------------
# 计算每个月的 Cycle 总数和离群 Cycle 的数量，用于绘制堆叠图
df_monthly_metrics <- df_metrics_extended %>%
  group_by(Year, Month) %>%
  summarise(
    # 确定该月是否是滑动窗口的初始阶段 (即 WL 统计量是否有 NA)
    Is_WL_Start = any(is.na(WL_D_SD_Median)),
    
    Total_Cycles = n_distinct(group.id),
    
    # dD 离群值计数和百分比
    Outlier_Cycles_D = sum(Outlier_D, na.rm = TRUE),
    Perc_Outlier_D = (Outlier_Cycles_D / Total_Cycles) * 100,
    
    # dO18 离群值计数和百分比
    Outlier_Cycles_O18 = sum(Outlier_O18, na.rm = TRUE),
    Perc_Outlier_O18 = (Outlier_Cycles_O18 / Total_Cycles) * 100,
    
    # 月度 Z-score 均值 (用于颜色映射)
    Mean_Z_D = mean(Z_Composite_D, na.rm = TRUE),
    Mean_Z_O18 = mean(Z_Composite_O18, na.rm = TRUE),
    
    # 获取该月任意 Cycle 的全局阈值 (用于透传)
    D_Threshold_SL = first(D_Threshold_SL),
    O18_Threshold_SL = first(O18_Threshold_SL),
    
    .groups = 'drop'
  ) %>%
  # 创建日期列用于绘图
  mutate(Date = as.Date(paste(Year, Month, "01", sep = "-")))

message("月度聚合计算完成。")

# ----------------------------------------------------------------------
# 6. 时间序列图绘制（Z-score 趋势）
# ----------------------------------------------------------------------

# 绘图函数：用于统一绘制 Z-score 时间序列图
plot_z_time_series <- function(data, z_col, threshold_col, y_label, title_label) {
  # 找出时间序列的最大日期用于放置文本标签
  Date_max <- max(data$Date, na.rm = TRUE)
  
  # 确保阈值是单个数值
  threshold_value <- unique(data[[threshold_col]])[1]
  
  # Z-score 的全局阈值 K_SL_SIGMA
  Z_threshold <- K_SL_SIGMA 

  p <- data %>%
    ggplot(aes(x = Date, y = !!sym(z_col))) +
    # 绘制 Z-score 趋势
    geom_line(color = "blue", alpha = 0.7) +
    geom_point(aes(size = 0.5), show.legend = FALSE) +
    # 标记离群值点 (Z-score > K_SL_SIGMA)
    geom_point(data = filter(data, !!sym(z_col) > Z_threshold), 
               aes(y = !!sym(z_col)), 
               color = "red", size = 1.5, shape = 21, fill = "red") +
    # 绘制阈值横线
    geom_hline(aes(yintercept = Z_threshold), 
               linetype = "dashed", color = "red", size = 0.5) +
    # 阈值横线添加文字注明
    geom_text(data = data.frame(Date = Date_max, 
                                Z_threshold = Z_threshold),
              aes(x = Date, y = Z_threshold, 
                  label = TeX(paste0("Z-score Threshold: ", K_SL_SIGMA, "$\\sigma$"))), 
              hjust = 1.05, vjust = -0.5, color = "red", size = 3) +
    
    labs(
      title = paste(title_label, "时间序列"),
      subtitle = paste("基于滑动窗口 (N=", WINDOW_WIDTH, ", s=", WINDOW_STEP, ") 和 全局 Z-score ($\\text{K}_{\\text{SL}}=$", K_SL_SIGMA, ") 综合判断"),
      x = "日期",
      y = y_label
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
  
  return(p)
}

# --- 绘制 Z_Composite_D 时间序列 ---
plot_Z_D <- plot_z_time_series(df_metrics_extended, "Z_Composite_D", "Z_Threshold_SL", 
                                 TeX("Z-Score ($\\text{dD}$)"), TeX("Z-Score of $\\text{dD}$ Composite Variability"))
ggsave(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_D.pdf"), plot_Z_D, width = 12, height = 6)
message("WVISS_Z_Composite_D.pdf 绘制完成。")

# --- 绘制 Z_Composite_O18 时间序列 (新增) ---
plot_Z_O18 <- plot_z_time_series(df_metrics_extended, "Z_Composite_O18", "Z_Threshold_SL", 
                                 TeX("Z-Score ($\\text{d}^{18}\\text{O}$)"), TeX("Z-Score of $\\text{d}^{18}\\text{O}$ Composite Variability"))
ggsave(file.path(VISUALIZATION_DIR, "WVISS_Z_Composite_O18.pdf"), plot_Z_O18, width = 12, height = 6)
message("WVISS_Z_Composite_O18.pdf 绘制完成。")


# ----------------------------------------------------------------------
# 7. 月度状态堆叠图绘制（离群值百分比）
# ----------------------------------------------------------------------

# 堆叠图函数：用于统一绘制月度离群值堆叠图
plot_stacked_discrete <- function(data, outlier_col, y_label, threshold_value) {
  
  # 用于判断离群值的列是 Outlier_D, Outlier_O18 或 Z_Composite_Combined (需要先计算)
  if (outlier_col == "Z_Composite_D") {
    perc_col <- "Perc_Outlier_D"
    title_var <- TeX("$\\text{dD}$")
  } else if (outlier_col == "Z_Composite_O18") {
    perc_col <- "Perc_Outlier_O18"
    title_var <- TeX("$\\text{d}^{18}\\text{O}$")
  } else if (outlier_col == "Z_Composite_Combined") {
    perc_col <- "Perc_Outlier_Combined" # 将在下面计算
    title_var <- TeX("Combined $\\text{dD}$ & $\\text{d}^{18}\\text{O}$")
  } else {
    stop("无效的 outlier_col 参数")
  }

  # 只有在 Combined 模式下需要计算新的百分比列
  if (outlier_col == "Z_Composite_Combined") {
    # 组合离群值：任一 Z-score 超过其对应的全局阈值，则标记为 1 (Outlier)
    data <- data %>%
      # 检查 Cycle 级别数据中是否有 D 或 O18 离群
      mutate(
        Z_Composite_Combined = as.integer(
          (Mean_Z_D > K_SL_SIGMA | Mean_Z_O18 > K_SL_SIGMA) & Total_Cycles > 0
        )
      ) %>%
      # 重新计算组合离群值的百分比
      group_by(Year, Month) %>%
      mutate(
        Outlier_Cycles_Combined = sum(df_metrics_extended %>% 
                                        filter(Year == first(Year), Month == first(Month)) %>%
                                        pull(Outlier_D) | df_metrics_extended %>% 
                                        filter(Year == first(Year), Month == first(Month)) %>%
                                        pull(Outlier_O18), na.rm = TRUE),
        Perc_Outlier_Combined = (Outlier_Cycles_Combined / Total_Cycles) * 100
      ) %>% ungroup()
  }
  
  # 确定用于颜色映射的 Z-score 列
  mean_z_col <- ifelse(outlier_col == "Z_Composite_D", "Mean_Z_D", "Mean_Z_O18")
  if (outlier_col == "Z_Composite_Combined") {
    # Combined 模式下，使用两个 Z-score 的最大值作为颜色映射
    data <- data %>% mutate(Mean_Z_Combined = pmax(Mean_Z_D, Mean_Z_O18, na.rm = TRUE))
    mean_z_col <- "Mean_Z_Combined"
  }


  p <- data %>%
    ggplot(aes(x = Date, y = !!sym(perc_col))) +
    # 堆叠条形图
    geom_bar(stat = "identity", position = "identity", fill = "gray90", color = "gray50", linewidth = 0.2) +
    # 叠加离群值百分比，颜色映射到该月的 Z-score 均值
    geom_bar(data = filter(data, !!sym(perc_col) > 0), 
             stat = "identity", position = "identity", 
             aes(fill = !!sym(mean_z_col)), color = "black", linewidth = 0.2) +
    
    # 颜色标尺
    scale_fill_gradientn(
      colors = c("green", "yellow", "orange", "red"),
      values = scales::rescale(c(0, K_SL_SIGMA * 0.5, K_SL_SIGMA * 0.8, K_SL_SIGMA)),
      limits = c(0, max(data[[mean_z_col]], na.rm = TRUE)),
      name = "Mean Cycle Z-score"
    ) +
    
    labs(
      title = paste(title_var, "月度周期状态堆叠图 (离群值百分比)"),
      subtitle = paste("阈值 Z-score >", K_SL_SIGMA, " (Median + 5*MAD)"),
      x = "日期",
      y = "离群周期百分比 (%)"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}

# --- 1. 绘制基于 Z_Composite_D 的堆叠图 ---
plot_stacked_D <- plot_stacked_discrete(df_monthly_metrics, "Z_Composite_D", 
                                        "$\\text{dD}$ Variability Status", df_global_stats$D_SD_Median)
ggsave(file.path(VISUALIZATION_DIR, "WVISS_Monthly_Status_Stacked_Discrete_D.pdf"), plot_stacked_D, width = 14, height = 7)
message("WVISS_Monthly_Status_Stacked_Discrete_D.pdf 绘制完成。")

# --- 2. 绘制基于 Z_Composite_O18 的堆叠图 (新增) ---
plot_stacked_O18 <- plot_stacked_discrete(df_monthly_metrics, "Z_Composite_O18", 
                                          "$\\text{d}^{18}\\text{O}$ Variability Status", df_global_stats$O18_SD_Median)
ggsave(file.path(VISUALIZATION_DIR, "WVISS_Monthly_Status_Stacked_Discrete_O18.pdf"), plot_stacked_O18, width = 14, height = 7)
message("WVISS_Monthly_Status_Stacked_Discrete_O18.pdf 绘制完成。")

# --- 3. 绘制基于 Z_Composite_Combined (D 或 O18 离群) 的堆叠图 (新增) ---
plot_stacked_Combined <- plot_stacked_discrete(df_monthly_metrics, "Z_Composite_Combined", 
                                               "Combined D & $\\text{d}^{18}\\text{O}$ Variability Status", df_global_stats$D_SD_Median)
ggsave(file.path(VISUALIZATION_DIR, "WVISS_Monthly_Status_Stacked_Discrete_Combined.pdf"), plot_stacked_Combined, width = 14, height = 7)
message("WVISS_Monthly_Status_Stacked_Discrete_Combined.pdf 绘制完成。")

message("\n所有分析和可视化文件已保存到：")
message(OUTPUT_DIR)
# ----------------------------------------------------------------------
# 8. 结果输出
# ----------------------------------------------------------------------
write.csv(df_metrics_extended, 
          file.path(OUTPUT_DIR, "WVISS_Cycle_Metrics_Extended.csv"), 
          row.names = FALSE)
write.csv(df_monthly_metrics, 
          file.path(OUTPUT_DIR, "WVISS_Monthly_Metrics.csv"), 
          row.names = FALSE)

message("最终结果 CSV 文件已保存。脚本运行结束。")