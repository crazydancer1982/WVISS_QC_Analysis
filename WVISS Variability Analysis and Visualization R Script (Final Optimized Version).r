# ======================================================================
# WVISS Variability Analysis and Visualization R Script (Section 1 & 2 Update)
# ======================================================================

# ----------------------------------------------------------------------
# 1. Load Required Libraries & QC 参数定义
# ----------------------------------------------------------------------
# Ensure these packages are installed: install.packages(c("dplyr", "tidyr", "zoo", "ggplot2", "scales", "RcppRoll"))
library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
library(scales) 
tryCatch(library(RcppRoll), error = function(e) message("RcppRoll not found. Using zoo::rollapply (slower) for rolling median."))
# 如果您在后续的完整脚本中使用了 maditr，请在此处也加载
# library(maditr) 

options("digits" = 9) 

# --- Global QC Parameters Setting (与完整脚本保持一致) ---
# 滑动窗口参数 (WL)
WINDOW_WIDTH <- 270L  # 滑动窗口宽度 N (Cycles) (更新为 270)
WINDOW_STEP <- 5L     # 滑动窗口步长 s (Cycles) (新增为 5)

# 全局（SL - Static/System Level）和 局部（WL - Window Level）离群值阈值
K_SL_SIGMA <- 5.0     # 全局离群值阈值的 Sigma 倍数 (K_MAD_THRESHOLD 更新为 5.0)
K_WL_SIGMA <- 3.0     # 局部离群值阈值的 Sigma 倍数 (新增/保持 3.0)

ISOTOPES <- c("D", "O18") # Isotopes

message(paste("QC 参数设置：K_SL_SIGMA =", K_SL_SIGMA, 
              ", K_WL_SIGMA =", K_WL_SIGMA, 
              ", N (Window Width) =", WINDOW_WIDTH, 
              ", s (Step) =", WINDOW_STEP))

# ----------------------------------------------------------------------
# 2. Data Loading and Pre-processing 
# ----------------------------------------------------------------------

# **ACTION REQUIRED: Manually set the path where WVISS.all.csv is located**
wd.synthesis <- '/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'
data_path <- file.path(wd.synthesis, 'WVISS.all.csv')

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

# 更改工作目录到 QC scripts 所在路径 (后续步骤可能需要)
wd.QC_A <- '/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/'
setwd(wd.QC_A)

# ----------------------------------------------------------------------
# 4. Calculation of Robust Metrics (Median, MAD, CL) and Outlier Flag
# ----------------------------------------------------------------------


# Ensure df_cycle_metrics is sorted for rolling window calculation
df_cycle_metrics <- df_cycle_metrics %>% arrange(Cycle_Seq) 

for (iso in ISOTOPES) {
  sd_col <- paste0(iso, "_SD")
  range_col <- paste0(iso, "_Range")
  
  # --- 4.1. Calculate 10-Year Robust Statistics (Median, MAD) ---
  median_sd <- median(df_cycle_metrics[[sd_col]], na.rm = TRUE)
  mad_sd <- mad(df_cycle_metrics[[sd_col]], constant = 1.4826, na.rm = TRUE)
  median_range <- median(df_cycle_metrics[[range_col]], na.rm = TRUE)
  mad_range <- mad(df_cycle_metrics[[range_col]], constant = 1.4826, na.rm = TRUE)
  
  # --- 4.2. Add Robust CL Thresholds (Median + K*MAD) and 10Y Median ---
  df_cycle_metrics <- df_cycle_metrics %>%
    mutate(
      # Baseline
      !!paste0(sd_col, "_Median10Y") := median_sd, 
      !!paste0(range_col, "_Median10Y") := median_range, 
      # Robust Control Limit (CL)
      !!paste0(sd_col, "_CL3MAD") := median_sd + K_MAD_THRESHOLD * mad_sd, 
      !!paste0(range_col, "_CL3MAD") := median_range + K_MAD_THRESHOLD * mad_range
    )
  
  # --- 4.3. Calculate 210 Cycle Rolling Median Trend ---
  if (exists("roll_median", where = as.environment("package:RcppRoll"))) {
      df_cycle_metrics <- df_cycle_metrics %>%
        mutate(
          !!paste0(sd_col, "_Trend") := RcppRoll::roll_median(.data[[sd_col]], n = ROLLING_WINDOW_N, fill = NA, align = "center", na.rm = TRUE),
          !!paste0(range_col, "_Trend") := RcppRoll::roll_median(.data[[range_col]], n = ROLLING_WINDOW_N, fill = NA, align = "center", na.rm = TRUE)
        )
  } else {
      df_cycle_metrics <- df_cycle_metrics %>%
        mutate(
          !!paste0(sd_col, "_Trend") := zoo::rollapply(.data[[sd_col]], width = ROLLING_WINDOW_N, FUN = median, by = 1, fill = NA, align = "center", na.rm = TRUE),
          !!paste0(range_col, "_Trend") := zoo::rollapply(.data[[range_col]], width = ROLLING_WINDOW_N, FUN = median, by = 1, fill = NA, align = "center", na.rm = TRUE)
        )
  }
  
  # --- 4.4. Calculate Z_Composite (Unified Standardized Volatility Index) ---
  z_sd_col <- paste0("Z_SD_", iso)
  z_range_col <- paste0("Z_Range_", iso)
  
  df_cycle_metrics <- df_cycle_metrics %>%
    mutate(
      !!z_sd_col := (.data[[sd_col]] - median_sd) / mad_sd,
      !!z_range_col := (.data[[range_col]] - median_range) / mad_range
    )
}

# Final Z_Composite: Max absolute value of the two Z-scores
df_cycle_metrics <- df_cycle_metrics %>%
  mutate(
    Z_Composite_D = abs(Z_SD_D) %>% pmax(abs(Z_Range_D)),
    Z_Composite_O18 = abs(Z_SD_O18) %>% pmax(abs(Z_Range_O18)),
    # --- 4.5. Outlier Flag Calculation ---
    D_Outlier = (D_SD > D_SD_CL3MAD) | (D_Range > D_Range_CL3MAD),
    O18_Outlier = (O18_SD > O18_SD_CL3MAD) | (O18_Range > O18_Range_CL3MAD),
    # Any_Outlier is TRUE if D or O18 is an outlier
    Any_Outlier = D_Outlier | O18_Outlier,
    # Handle NA Outlier flag by setting to FALSE
    Any_Outlier = replace_na(Any_Outlier, FALSE)
  )

message("Robust metrics and Z_Composite calculation completed.")


# ----------------------------------------------------------------------
# 5. Monthly Aggregation: Calculate Exceedance Rate (MUST BE EXECUTED)
# ----------------------------------------------------------------------

df_monthly_status <- df_cycle_metrics %>%
  # Filter out cycles with missing Year/Month information
  filter(!is.na(Year) & !is.na(Month)) %>% 
  
  # Aggregate by Year and Month
  group_by(Year, Month) %>%
  summarise(
    Total_Cycles = n(),
    # Calculate Percentage of Outliers
    Perc_Outlier = sum(Any_Outlier) / Total_Cycles * 100,
    .groups = 'drop'
  ) %>%
  
  # Create Date and Discrete X-axis Label (YYYY-MM)
  mutate(
    Date = as.Date(paste(Year, Month, "01", sep = "-")),
    Year_Month_Label = format(Date, "%Y-%m"),
    Perc_Outlier = replace_na(Perc_Outlier, 0)
  ) %>%
  arrange(Date)

# Calculate yearly status for the subtitle (Annual Average)
df_yearly_status <- df_cycle_metrics %>%
  filter(!is.na(Year)) %>%
  group_by(Year) %>%
  summarise(
    Total_Cycles = n(),
    Perc_Outlier = sum(Any_Outlier) / Total_Cycles * 100,
    .groups = 'drop'
  )

message(paste("Monthly aggregation table df_monthly_status created. Total monthly data points:", nrow(df_monthly_status)))


# ----------------------------------------------------------------------
# 6. Visualization: Output to Separate PDF Files (English Labels)
# ----------------------------------------------------------------------

k_threshold <- K_MAD_THRESHOLD
roll_n <- ROLLING_WINDOW_N

message("\nGenerating individual PDF files with English labels...")

# --- 6.1 & 6.2. SD and Range Time Series Plots (Loop for D and O18) ---

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
  
  # Save SD Plot
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
  
  # Save Range Plot
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

# Save Z_Composite Plot
pdf_filename_zcomp <- "WVISS_Z_Composite_D.pdf"
pdf(pdf_filename_zcomp, width = 11, height = 8.5); print(p_zcomp); dev.off()
message(paste("- Saved:", pdf_filename_zcomp))


# ----------------------------------------------------------------------
# 6.4. Monthly Stacked Column Chart (FINAL FIX: Discrete X-axis for 109 bars)
# ----------------------------------------------------------------------

# 1. Data Transformation for Stacking
df_monthly_stacked <- df_monthly_status %>%
    mutate(Perc_Non_Outlier = 100 - Perc_Outlier) %>%
    
    tidyr::pivot_longer(
      cols = c(Perc_Outlier, Perc_Non_Outlier),
      names_to = "Status_Type_Raw",
      values_to = "Percentage"
    ) %>%
    mutate(
      Status_Type = factor(Status_Type_Raw, 
                           levels = c("Perc_Non_Outlier", "Perc_Outlier"),
                           labels = c("Non-Exceedance", "Exceedance"))
    )

# 2. Define Colors and Subtitle
STACKED_COLORS <- c("Non-Exceedance" = "#009E73", "Exceedance" = "#D55E00")
annual_avg_exceedance <- round(mean(df_yearly_status$Perc_Outlier), 2)

# 3. Plotting
p_monthly_status_stacked <- df_monthly_stacked %>%
  # Use discrete X-axis (Year_Month_Label) to ensure all 109 months are rendered as separate bars
  ggplot(aes(x = Year_Month_Label, y = Percentage, fill = Status_Type)) + 
  
  geom_col(position = "stack") + 
  
  # Color Mapping
  scale_fill_manual(
    name = "Cycle Precision Status",
    values = STACKED_COLORS
  ) +
  
  # X-axis Optimization: Show only yearly labels
  scale_x_discrete(
    breaks = df_monthly_stacked$Year_Month_Label[grep("-01$", df_monthly_stacked$Year_Month_Label)]
  ) +
  
  # Labels and Theme
  labs(
    title = "Monthly Cycle Precision Status: Exceedance vs. Non-Exceedance (Discrete Axis)",
    subtitle = paste0("Annual Average Exceedance Rate: ", annual_avg_exceedance, "%"),
    y = "Cycle Percentage (%)",
    x = "Time (Year-Month)"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)), 
    limits = c(0, 100),
    breaks = c(0, 5, 10, 50, 100)
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        # Rotate X-axis labels to prevent overlap
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        # Add minor grid lines for better month alignment reference
        panel.grid.minor.x = element_line(color = "gray90", linetype = "dotted")) 

# Save Stacked Plot
pdf_filename_monthly_stacked_discrete <- "WVISS_Monthly_Status_Stacked_Discrete.pdf"
pdf(pdf_filename_monthly_stacked_discrete, width = 11, height = 8.5); print(p_monthly_status_stacked); dev.off()
message(paste("- Saved:", pdf_filename_monthly_stacked_discrete))

message("\nAll 6 plots have been successfully generated into separate PDF files in the working directory.")
