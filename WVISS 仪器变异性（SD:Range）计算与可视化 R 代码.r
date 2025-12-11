# ======================================================================
# WVISS 仪器变异性（SD/Range）计算与可视化 R 代码 (增强版 - PDF 输出)
# 目标: 计算 Cycle-level SD/Range, 10年鲁棒阈值, Z_Composite, 月/年超限比例
# ======================================================================

# 1. 加载所需的库
# 确保已安装这些包：install.packages(c("dplyr", "tidyr", "zoo", "ggplot2", "scales"))
#install.packages("RcppRoll")
library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
library(scales) 
library(RcppRoll)
options("digits" = 9) 

# --- 全局参数设定 ---
ROLLING_WINDOW_N <- 210 # 滑动窗口大小 (Cycles)
ISOTOPES <- c("D", "O18") # 同位素种类
DELTA_COLS <- c("Ddel_mean", "O18del_mean") 
names(DELTA_COLS) <- ISOTOPES
# 鲁棒 CL 阈值的 K 值 (例如 3*MAD 对应 K=3)
K_MAD_THRESHOLD <- 3 

# ----------------------------------------------------------------------
# 2. 实际数据加载与预处理
# ----------------------------------------------------------------------
wd.synthesis <- '/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'
setwd(wd.synthesis)

# **请手动将此路径修改为 WVISS.all.csv 所在的实际目录**
data_path <- 'WVISS.all.csv' 
WVISS.all <- NULL
if (file.exists(data_path)) {
  WVISS.all <- read.csv(data_path, header = T, stringsAsFactors = FALSE)
  WVISS.all <- as_tibble(WVISS.all) 
  message("WVISS.all.csv 数据加载成功。")
} else {
  stop(paste("错误：WVISS.all.csv 文件不存在于路径：", data_path, ". 请检查路径和文件名。\n"))
}

wd.QC_A <- '/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/'
setwd(wd.QC_A)

# --- 关键修正：基于 group_DOY_mean 排序并定义 Cycle_Seq ---
WVISS.data <- WVISS.all %>%
  # 1. 计算每个 Cycle 的平均时间 (group_DOY_mean)
  group_by(group.id) %>%
  mutate(
    group_DOY_mean = mean(DOY_mean, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  # 2. 按 Cycle 平均时间排序数据框，以确定正确的 group.id 顺序
  arrange(group_DOY_mean)

# 3. 提取排序后的 group.id 列表，作为因子的 levels
chronological_groups <- unique(WVISS.data$group.id)
head(chronological_groups)

# 4. 创建最终的 WVISS.data 数据框，并生成 Cycle_Seq
WVISS.data <- WVISS.data %>%
  mutate(
    # 基于 2011-01-01 和 DOY_mean (连续天数) 计算 Date
    Start_Date = as.Date("2011-01-01"),
    Date = Start_Date + floor(DOY_mean) - 1,
    # *** 修正 Cycle 序号生成 ***：使用排序后的列表作为因子 levels
    Cycle_Seq = as.numeric(factor(group.id, levels = chronological_groups)),
    # 提取 Month 和 Year 用于分组分析
    Month = as.numeric(format(Date, "%m")),
    Year = as.numeric(format(Date, "%Y"))
  ) %>%
  select(-Start_Date)
str(WVISS.data)
head(WVISS.data$Date)
# ----------------------------------------------------------------------
# 3. Cycle 级别变异性（SD 和 Range）计算
# (保留多列分组以携带元数据到下一步)
# ----------------------------------------------------------------------
?group_by
df_cycle_metrics <- WVISS.data %>%
  # 保持多列分组：确保 Date, Year, Month, Cycle_Seq 都被保留
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
  # *** 关键修正：显式地按 Cycle_Seq 排序 ***
  arrange(Cycle_Seq) %>%
  # 移除 SD 或 Range 为 Inf 的行
  filter_at(vars(D_SD, D_Range, O18_SD, O18_Range), all_vars(!is.infinite(.)))

message(paste("Cycle 级别 SD/Range 计算完成。共计", nrow(df_cycle_metrics), "个 Cycle。"))

View(df_cycle_metrics)
dim(df_cycle_metrics)
dim(WVISS.all)/5

# ----------------------------------------------------------------------
# 4. 变异性指标计算：鲁棒基线、MAD 和 Z_Composite
# ----------------------------------------------------------------------

# 存储所有鲁棒统计量，用于 Z_Composite 计算
robust_stats <- list()

for (iso in ISOTOPES) {
  sd_col <- paste0(iso, "_SD")
  range_col <- paste0(iso, "_Range")
  
  # --- 4.1. 计算 10 年鲁棒统计量 (Median, MAD) ---
  median_sd <- median(df_cycle_metrics[[sd_col]], na.rm = TRUE)
  mad_sd <- mad(df_cycle_metrics[[sd_col]], constant = 1.4826, na.rm = TRUE)
  median_range <- median(df_cycle_metrics[[range_col]], na.rm = TRUE)
  mad_range <- mad(df_cycle_metrics[[range_col]], constant = 1.4826, na.rm = TRUE)
  
  # 存储统计量
  robust_stats[[iso]] <- list(
    SD_Median = median_sd,
    SD_MAD = mad_sd,
    Range_Median = median_range,
    Range_MAD = mad_range
  )
  
  # --- 4.2. 添加鲁棒 CL 阈值 (Median + 3*MAD) 和 Median10Y ---
  df_cycle_metrics <- df_cycle_metrics %>%
    mutate(
      !!paste0(sd_col, "_Median10Y") := median_sd, # 10年 Median (基线)
      !!paste0(sd_col, "_CL3MAD") := median_sd + K_MAD_THRESHOLD * mad_sd, # Median + 3*MAD (鲁棒 CL)
      !!paste0(range_col, "_Median10Y") := median_range, 
      !!paste0(range_col, "_CL3MAD") := median_range + K_MAD_THRESHOLD * mad_range
    )
  
  # --- 4.3. 计算 210 Cycle 滑动中位数趋势 (动态 CL 趋势) ---
  df_cycle_metrics <- df_cycle_metrics %>% arrange(Cycle_Seq) %>%
  mutate(
    # 使用 RcppRoll 替代 zoo::rollapply
    !!paste0(sd_col, "_Trend") := RcppRoll::roll_median(
      .data[[sd_col]], n = ROLLING_WINDOW_N, fill = NA, align = "center", na.rm = TRUE
    ),
    !!paste0(range_col, "_Trend") := RcppRoll::roll_median(
      .data[[range_col]], n = ROLLING_WINDOW_N, fill = NA, align = "center", na.rm = TRUE
    )
  )
  
  # --- 4.4. 计算 Z_Composite (统一标准化波动指标) ---
  # Z_SD 和 Z_Range
  z_sd_col <- paste0("Z_SD_", iso)
  z_range_col <- paste0("Z_Range_", iso)
  
  df_cycle_metrics <- df_cycle_metrics %>%
    mutate(
      !!z_sd_col := (.data[[sd_col]] - median_sd) / mad_sd,
      !!z_range_col := (.data[[range_col]] - median_range) / mad_range
    )
}

# 最终 Z_Composite: 两个 Z-score 中绝对值的最大值
df_cycle_metrics <- df_cycle_metrics %>%
  mutate(
    Z_Composite_D = abs(Z_SD_D) %>% pmax(abs(Z_Range_D)),
    Z_Composite_O18 = abs(Z_SD_O18) %>% pmax(abs(Z_Range_O18))
  )
View(df_cycle_metrics)
dim(df_cycle_metrics)
head(df_cycle_metrics)
message("鲁棒统计量和 Z_Composite 计算完成。")

# ----------------------------------------------------------------------
# 5. 月度/年度 超限百分比分析
# ----------------------------------------------------------------------

# 定义超限判定列：SD 或 Range 超过 Median + 3*MAD 阈值
df_cycle_metrics <- df_cycle_metrics %>%
  mutate(
    D_Outlier = (D_SD > D_SD_CL3MAD) | (D_Range > D_Range_CL3MAD),
    O18_Outlier = (O18_SD > O18_SD_CL3MAD) | (O18_Range > O18_Range_CL3MAD),
    Any_Outlier = D_Outlier | O18_Outlier
  )

# 计算月度和年度超限百分比
df_monthly_status <- df_cycle_metrics %>%
  group_by(Year, Month) %>%
  summarise(
    Total_Cycles = n(),
    Perc_Outlier = sum(Any_Outlier) / Total_Cycles * 100,
    .groups = 'drop'
  ) %>%
  mutate(Date = as.Date(paste(Year, Month, "01", sep = "-")))

df_yearly_status <- df_cycle_metrics %>%
  group_by(Year) %>%
  summarise(
    Total_Cycles = n(),
    Perc_Outlier = sum(Any_Outlier) / Total_Cycles * 100,
    .groups = 'drop'
  )

message("月度/年度超限百分比计算完成。")

# ----------------------------------------------------------------------
# 6. 可视化：输出为多文件 PDF (每图一文件，标签英文化)
# ----------------------------------------------------------------------

# 确保 K_MAD_THRESHOLD 和 ROLLING_WINDOW_N 已在全局设置中定义
k_threshold <- K_MAD_THRESHOLD
roll_n <- ROLLING_WINDOW_N

message("\nGenerating individual PDF files with English labels...")

# --- 6.1 & 6.2. SD 和 Range 时间序列图 (在 Loop 中生成) ---

for (iso in ISOTOPES) {
  sd_col <- paste0(iso, "_SD")
  range_col <- paste0(iso, "_Range")
  median_sd_col <- paste0(sd_col, "_Median10Y")
  cl3mad_sd_col <- paste0(sd_col, "_CL3MAD")
  trend_sd_col <- paste0(sd_col, "_Trend")
  
  # --- 6.1. SD Time Series Plot (p_sd) ---
  p_sd <- df_cycle_metrics %>%
    ggplot(aes(x = Date)) + 
    geom_point(aes(y = .data[[sd_col]], 
                   color = .data[[sd_col]] > .data[[cl3mad_sd_col]]), 
               size = 0.8, alpha = 0.5) +
    geom_line(aes(y = .data[[median_sd_col]], color = "Median (Baseline)"), 
              linetype = "solid", size = 0.5) +
    geom_line(aes(y = .data[[cl3mad_sd_col]], color = "Median + 3*MAD (Robust CL)"), 
              linetype = "dashed", size = 0.8) +
    geom_line(aes(y = .data[[trend_sd_col]], color = "Rolling Median Trend"), 
              size = 1.0) +
    scale_color_manual(
      name = "Variability Status and Trend",
      values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2", "Median (Baseline)" = "#009E73", 
                 "Median + 3*MAD (Robust CL)" = "#E69F00", "Rolling Median Trend" = "#CC79A7"),
      labels = c("TRUE" = paste("SD > Median +", k_threshold, "*MAD"), 
                 "FALSE" = "SD ≤ Threshold",
                 "Median (Baseline)" = "10-Year Absolute Median (Baseline)",
                 "Median + 3*MAD (Robust CL)" = paste("Robust Control Limit (3*MAD)"),
                 "Rolling Median Trend" = paste0(roll_n, " Cycle Rolling Median (Dynamic Trend)")
    )
    ) +
    labs(
      title = paste(iso, "SD Variability Time Series: Absolute Fluctuation and Dynamic Trend"),
      y = paste(iso, "SD (Cycle-level)"),
      x = "Date"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # 输出 SD 图到文件
  pdf_filename_sd <- paste0("WVISS_", iso, "_SD_TimeSeries.pdf")
  pdf(pdf_filename_sd, width = 11, height = 8.5); print(p_sd); dev.off()
  message(paste("- Saved:", pdf_filename_sd))
  
  # --- 6.2. Range Time Series Plot (p_range) ---
  p_range <- df_cycle_metrics %>%
    ggplot(aes(x = Date)) +
    geom_point(aes(y = .data[[range_col]], 
                   color = .data[[range_col]] > .data[[paste0(range_col, "_CL3MAD")]]), 
               size = 0.8, alpha = 0.5) +
    geom_line(aes(y = .data[[paste0(range_col, "_Median10Y")]], color = "Median (Baseline)"), 
              linetype = "solid", size = 0.5) +
    geom_line(aes(y = .data[[paste0(range_col, "_CL3MAD")]], color = "Median + 3*MAD (Robust CL)"), 
              linetype = "dashed", size = 0.8) +
    geom_line(aes(y = .data[[paste0(range_col, "_Trend")]], color = "Rolling Median Trend"), 
              size = 1.0) +
    scale_color_manual(
      name = "Variability Status and Trend",
      values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2", "Median (Baseline)" = "#009E73", 
                 "Median + 3*MAD (Robust CL)" = "#E69F00", "Rolling Median Trend" = "#CC79A7"),
      labels = c("TRUE" = paste("Range > Median +", k_threshold, "*MAD"), 
                 "FALSE" = "Range ≤ Threshold",
                 "Median (Baseline)" = "10-Year Absolute Median (Baseline)",
                 "Median + 3*MAD (Robust CL)" = paste("Robust Control Limit (3*MAD)"),
                 "Rolling Median Trend" = paste0(roll_n, " Cycle Rolling Median (Dynamic Trend)")
    )
    ) +
    labs(
      title = paste(iso, "Range Variability Time Series: More Sensitive to Extremes"),
      y = paste(iso, "Range (Cycle-level)"),
      x = "Date"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # 输出 Range 图到文件
  pdf_filename_range <- paste0("WVISS_", iso, "_Range_TimeSeries.pdf")
  pdf(pdf_filename_range, width = 11, height = 8.5); print(p_range); dev.off()
  message(paste("- Saved:", pdf_filename_range))
}


# --- 6.3. Z_Composite Control Chart (dD) ---
p_zcomp <- df_cycle_metrics %>%
  ggplot(aes(x = Date, y = Z_Composite_D)) +
  geom_point(aes(color = Z_Composite_D > k_threshold), size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = k_threshold, linetype = "dashed", color = "#D55E00", size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 0.5) +
  scale_color_manual(
    name = "Outlier Status",
    values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
    labels = c("TRUE" = paste("Z >", k_threshold, "(Cycle Precision Out of Control)"), "FALSE" = paste("Z ≤", k_threshold))
  ) +
  labs(
    title = "dD Z-Composite Robust Control Chart",
    subtitle = paste("Z-Composite is the max absolute value of robustly standardized SD and Range. Threshold K=", k_threshold, "."),
    y = "Z_Composite (dD)",
    x = "Date"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# 输出 Z_Composite 图到文件
pdf_filename_zcomp <- "WVISS_Z_Composite_D.pdf"
pdf(pdf_filename_zcomp, width = 11, height = 8.5); print(p_zcomp); dev.off()
message(paste("- Saved:", pdf_filename_zcomp))
# --- 6.4. Monthly Exceedance Rate Chart (FINAL: Stacked Column Chart) ---

# 1. 数据重塑：计算 Non-Exceedance 比例并转换为长格式 (Tidy Data)
df_monthly_stacked <- df_monthly_status %>%
    # 增加 Non-Exceedance 比例列
    mutate(Perc_Non_Outlier = 100 - Perc_Outlier) %>%
    
    # 转换为长格式：将 Exceedance 和 Non-Exceedance 比例转为两行
    tidyr::pivot_longer(
      cols = c(Perc_Outlier, Perc_Non_Outlier),
      names_to = "Status_Type_Raw",
      values_to = "Percentage"
    ) %>%
    
    # 定义因子顺序和标签，确保 Non-Exceedance 在底部 (绿色)
    mutate(
      Status_Type = factor(Status_Type_Raw, 
                           levels = c("Perc_Non_Outlier", "Perc_Outlier"),
                           labels = c("Non-Exceedance", "Exceedance"))
    )

# 2. 定义颜色（简化为两色：绿色/红色）
STACKED_COLORS <- c("Non-Exceedance" = "#009E73", "Exceedance" = "#D55E00")

# 提取年度平均超限比例作为副标题
annual_avg_exceedance <- round(mean(df_yearly_status$Perc_Outlier), 2)

# 3. 绘制堆叠柱状图
p_monthly_status_stacked <- df_monthly_stacked %>%
  # 填充颜色 (fill) 必须与 Status_Type 绑定，position="stack" 确保堆叠
  ggplot(aes(x = Date, y = Percentage, fill = Status_Type)) + 
  
  geom_col(position = "stack", width = 28) + 
  
  # 颜色映射
  scale_fill_manual(
    name = "Cycle Precision Status",
    values = STACKED_COLORS
  ) +
  
  # 2. Management Threshold Lines (这里不再是 CL 阈值线，而是 Exceedance Rate 的管理阈值线)
  # 注意: 在堆叠图上绘制阈值线 (5%, 10%) 需要特殊的计算，但由于图表目标是展示 100%，
  #       我们只保留 100% 顶线，并移除 5%/10% 的水平线，改为强调颜色编码。
  
  # 确保 Y 轴上限为 100%
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)), 
    limits = c(0, 100),
    breaks = c(0, 5, 10, 50, 100) # 强调 5% 和 10% 的管理阈值
  ) +
  
  # X轴优化: 增加次刻度线
  scale_x_date(
    date_breaks = "1 year", 
    date_labels = "%Y",
    date_minor_breaks = "3 month"
  ) +
  
  # 5. Labels and Theme
  labs(
    title = "Monthly Cycle Precision Status: Exceedance vs. Non-Exceedance",
    subtitle = paste0("Annual Average Exceedance Rate: ", annual_avg_exceedance, "%"),
    y = "Cycle Percentage (%)",
    x = "Time (Monthly Aggregation Date)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid.minor.x = element_line(color = "gray90", linetype = "dotted"))

# 输出堆叠图到文件
pdf_filename_monthly_stacked <- "WVISS_Monthly_Status_Stacked.pdf"
pdf(pdf_filename_monthly_stacked, width = 11, height = 8.5); print(p_monthly_status_stacked); dev.off()
message(paste("- Saved:", pdf_filename_monthly_stacked))

message("\nAll 6 plots have been successfully generated into separate PDF files in the working directory.")


# ----------------------------------------------------------------------
# 诊断：检查月度数据点密度
# ----------------------------------------------------------------------

# 1. 获取 df_monthly_status 中的最小和最大日期
start_date <- min(df_monthly_status$Date)
end_date <- max(df_monthly_status$Date)

# 2. 生成所有预期月份的序列 (从起始月到结束月，每月 1 号)
# 确保 zoo 包已加载 (在前面的步骤中已加载)
all_months <- seq.Date(from = start_date, to = end_date, by = "month")

# 3. 统计实际和预期的月份数量
actual_data_months <- nrow(df_monthly_status)
total_expected_months <- length(all_months)

message(paste("--- 月度数据连续性检查 ---"))
message(paste("数据起始日期:", start_date, "结束日期:", end_date))
message(paste("预期总月数 (包括空月):", total_expected_months))
message(paste("实际有数据的月数 (df_monthly_status行数):", actual_data_months))
message(paste("数据覆盖率:", round(actual_data_months / total_expected_months * 100, 2), "%"))

# 4. 检查缺失的月份
missing_months <- all_months[!all_months %in% df_monthly_status$Date]

if (length(missing_months) > 0) {
    message(paste("--- 缺失数据月份数量:", length(missing_months), "---"))
    print("前 10 个缺失数据月份 (Month of 01):")
    print(head(missing_months, 10))
} else {
    message("所有月份都有数据，稀疏性是视觉压缩或数据点数量过少的误解。")
}

# ----------------------------------------------------------------------
# 打印 df_monthly_status 的头部和尾部，确认 Year/Month 范围
# ----------------------------------------------------------------------
message("\n--- df_monthly_status (Head) ---")
print(head(df_monthly_status))
message("\n--- df_monthly_status (Tail) ---")
print(tail(df_monthly_status))
