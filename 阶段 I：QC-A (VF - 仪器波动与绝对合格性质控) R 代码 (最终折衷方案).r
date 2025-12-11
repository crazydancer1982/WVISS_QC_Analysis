# ----------------------------------------------------------------------
# 阶段 I：QC-A (VF - 仪器波动与绝对合格性质控) R 代码 (最终方案 N=270, k=7, SL K=5, CL K=3)
# **更新：SL K 阈值改为 5，以降低对 O18 的敏感度**
# ----------------------------------------------------------------------
wd.QC_A <- '/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/'
setwd(wd.QC_A)
# ----------------------------------------------------------------------
# 阶段 I：QC-A (VF - 仪器波动与绝对合格性质控) R 代码 (增强版)
# 目标: 优化 QC 统计分析，区分 CL/SL Outlier 的数量和 Cycle ID 重叠情况。
# 参数: N_WINDOW=270, K_CL_SIGMA=3, K_SL_SIGMA=5, K_STEP=7
# ----------------------------------------------------------------------

# 加载所需的库
# 确保已安装：dplyr, tidyr, zoo
library(dplyr)
library(tidyr)
library(zoo)

options("digits" = 9) 

# --- 全局参数设定 (最终确定方案) ---
N_WINDOW <- 270     # 窗口长度 N* (保持 270)
K_CL_SIGMA <- 3     # 动态 CL 阈值 K (保持 3)
K_SL_SIGMA <- 5     # 静态 SL 阈值 K (更新为 5)
K_STEP <- 5         # 步长 k* (保持 7)

# MAD 鲁棒性函数 (重新定义，确保可用)
robust_mad <- function(x) {
  # 确保在计算 median 和 mad 时处理 NA
  mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
}

# ----------------------------------------------------------------------
# 1. 实际数据加载与预处理 (沿用现有代码逻辑)
# ----------------------------------------------------------------------

# **请手动将此路径修改为 WVISS.all.csv 所在的实际目录**
wd.synthesis <- '/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'

WVISS.all <- NULL
if (dir.exists(wd.synthesis)) {
  setwd(wd.synthesis)
  if (file.exists('WVISS.all.csv')) {
    WVISS.all <- read.csv('WVISS.all.csv', header = T, stringsAsFactors = FALSE)
    WVISS.all <- as_tibble(WVISS.all) 
    message("WVISS.all.csv 数据加载成功。")
  } else {
    stop("错误：WVISS.all.csv 文件不存在于指定路径。请检查 wd.synthesis 和文件名。\n")
  }
} else {
  stop("错误：工作目录不存在。请在运行前手动设置 wd.synthesis 为正确路径！\n")
}
wd.QC_A <- '/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/'
setwd(wd.QC_A)

# --- 1.1 核心偏差计算 ---
RD_smow <- 0.00015576
RO_smow <- 0.0020052

WVISS.all <- WVISS.all %>%
  mutate(
    dO_obs = ((O18ppm_mean / O16ppm_mean) / RO_smow - 1) * 1000,
    dD_obs = ((HODppm_mean / (O16ppm_mean * 2)) / RD_smow - 1) * 1000,
    dO_raw_bias = dO_obs - o_wahaha,
    dD_raw_bias = dD_obs - d_wahaha,
    year = substr(group.id, 1, 4), 
    month = substr(group.id, 5, 6),
    year_month = substr(group.id, 1, 6)
  ) %>%
  select(group.id, Gradient, dD_raw_bias, dO_raw_bias, year, month, year_month) 

# --- 1.2 数据重塑：长格式转宽格式 ---
df_wide <- WVISS.all %>%
  pivot_wider(
    id_cols = c(group.id, year, month, year_month),
    names_from = Gradient,
    values_from = c(dD_raw_bias, dO_raw_bias),
    names_prefix = "G"
  ) %>%
  arrange(group.id) %>%
  mutate(Cycle_Index = row_number()) %>% 
  rename(Cycle_ID = group.id)

# ----------------------------------------------------------------------
# 2. 准备 QC-A 核心指标数据：SD 和 Range
# ----------------------------------------------------------------------

isotopes <- c("D", "O18")
G_cols <- list(
  D = paste0("dD_raw_bias_G", 1:5),
  O18 = paste0("dO_raw_bias_G", 1:5)
)

df_qc_a <- df_wide 

for (iso in isotopes) {
  df_qc_a[[paste0("SD_", iso)]] <- apply(df_qc_a[, G_cols[[iso]]], 1, sd, na.rm = TRUE)
  
  if (iso == "D") {
    df_qc_a[[paste0("Range_", iso)]] <- abs(df_qc_a[[paste0("dD_raw_bias_G", 5)]] - df_qc_a[[paste0("dD_raw_bias_G", 1)]])
  } else { 
    df_qc_a[[paste0("Range_", iso)]] <- abs(df_qc_a[[paste0("dO_raw_bias_G", 5)]] - df_qc_a[[paste0("dO_raw_bias_G", 1)]])
  }
}

# ----------------------------------------------------------------------
# 3. 计算 Z_Composite 与静态规格限 (SL)
# ** 注意：使用 K_SL_SIGMA = 5 **
# ----------------------------------------------------------------------

w1 <- 0.4 
w2 <- 0.6 

for (iso in isotopes) {
  SD_col <- paste0("SD_", iso)
  Range_col <- paste0("Range_", iso)
  
  # (a) 计算标准化 Z 分数的基线统计量 (基于全部数据)
  mu_sd_all <- median(df_qc_a[[SD_col]], na.rm = TRUE)
  sigma_sd_all <- robust_mad(df_qc_a[[SD_col]])
  mu_range_all <- median(df_qc_a[[Range_col]], na.rm = TRUE)
  sigma_range_all <- robust_mad(df_qc_a[[Range_col]])
  
  # (b) 计算标准化 Z 分数
  df_qc_a[[paste0("Z_SD_", iso)]] <- (df_qc_a[[SD_col]] - mu_sd_all) / sigma_sd_all
  df_qc_a[[paste0("Z_Range_", iso)]] <- (df_qc_a[[Range_col]] - mu_range_all) / sigma_range_all
  
  # (c) 计算 Z_Composite
  Z_Comp_col <- paste0("Z_Composite_", iso)
  df_qc_a[[Z_Comp_col]] <- w1 * df_qc_a[[paste0("Z_SD_", iso)]] + w2 * df_qc_a[[paste0("Z_Range_", iso)]]
  
  # (d) 静态规格限 SL (基于全部 Z_Composite 数据)
  # SL = 全部数据 Z_Composite 的中位数 + K_SL_SIGMA * MAD
  mu_z_all <- median(df_qc_a[[Z_Comp_col]], na.rm = TRUE)
  sigma_z_all <- robust_mad(df_qc_a[[Z_Comp_col]])
  
  # *** 核心修改点：使用 K_SL_SIGMA (5) ***
  df_qc_a[[paste0("SL_", iso)]] <- mu_z_all + K_SL_SIGMA * sigma_z_all 
}

# ----------------------------------------------------------------------
# 4. 核心：计算动态控制限 (CL) 和最终 QC 判定
# **注意：使用 K_CL_SIGMA = 3**
# ----------------------------------------------------------------------

for (iso in isotopes) {
  Z_Comp_col <- paste0("Z_Composite_", iso)
  
  # --- 4.1 计算滚动统计量 (N=270 窗口) ---
  mu_roll <- rollapply(
    data = df_qc_a[[Z_Comp_col]], 
    width = N_WINDOW, 
    FUN = median, 
    align = "right", 
    fill = NA  
  )
  df_qc_a[[paste0("Mu_Roll_", iso)]] <- mu_roll

  sigma_roll <- rollapply(
    data = df_qc_a[[Z_Comp_col]], 
    width = N_WINDOW, 
    FUN = robust_mad, 
    align = "right", 
    fill = NA 
  )
  df_qc_a[[paste0("Sigma_Roll_", iso)]] <- sigma_roll

  # --- 4.2 计算动态控制限 (CL) ---
  CL_col <- paste0("CL_", iso)
  # Dynamic CL = Rolling Median + K_CL_SIGMA (3) * Rolling MAD 
  df_qc_a[[CL_col]] <- df_qc_a[[paste0("Mu_Roll_", iso)]] + K_CL_SIGMA * df_qc_a[[paste0("Sigma_Roll_", iso)]]
  
  # --- 4.3 最终 QC 判定 (CL 和 SL) ---
  
  # 1. 动态 CL Outlier (z_conposite_D/18O_CL) - 使用 K_CL_SIGMA (3) 和 K_STEP (7)
  CL_Outlier_col <- paste0("z_conposite_", iso, "_CL")
  df_qc_a[[CL_Outlier_col]] <- FALSE
  
  check_cycles <- which(df_qc_a$Cycle_Index %% K_STEP == 0)
  valid_check_cycles <- check_cycles[check_cycles >= N_WINDOW] 
  
  df_qc_a[[CL_Outlier_col]][valid_check_cycles] <- 
    df_qc_a[[Z_Comp_col]][valid_check_cycles] > df_qc_a[[CL_col]][valid_check_cycles]
  
  # 2. 全局 SL Outlier (z_conposite_D/18O_SL) - 使用 K_SL_SIGMA (5)
  SL_Outlier_col <- paste0("z_conposite_", iso, "_SL")
  df_qc_a[[SL_Outlier_col]] <- df_qc_a[[Z_Comp_col]] > df_qc_a[[paste0("SL_", iso)]]
}

# 3. 最终 Cycle 判定：任一情况满足，则标记为 Outlier
df_qc_a$outlier <- 
  df_qc_a$z_conposite_D_CL | 
  df_qc_a$z_conposite_D_SL | 
  df_qc_a$z_conposite_O18_CL | 
  df_qc_a$z_conposite_O18_SL

# ----------------------------------------------------------------------
# 5. 增强统计：CL 与 SL Outlier 的集合分析 (交集/差集)
# ----------------------------------------------------------------------

# (a) 定义总 CL 异常和总 SL 异常 (D 或 O18 任一超限)
df_qc_a <- df_qc_a %>%
  mutate(
    Total_CL_Outlier = z_conposite_D_CL | z_conposite_O18_CL,
    Total_SL_Outlier = z_conposite_D_SL | z_conposite_O18_SL
  )

# (b) 集合操作
df_qc_a <- df_qc_a %>%
  mutate(
    CL_AND_SL_Outlier = Total_CL_Outlier & Total_SL_Outlier,
    SL_ONLY_Outlier = Total_SL_Outlier & !Total_CL_Outlier,
    CL_ONLY_Outlier = Total_CL_Outlier & !Total_SL_Outlier
  )

# --- 5.1 按年份统计 (总 Outlier) ---
df_yearly_stats <- df_qc_a %>%
  group_by(year) %>%
  summarise(
    Total_Cycles = n(),
    Outlier_Count = sum(outlier, na.rm = TRUE),
    Normal_Count = Total_Cycles - Outlier_Count
  ) %>%
  mutate(
    Outlier_Rate = round(Outlier_Count / Total_Cycles * 100, 2)
  )

# --- 5.2 SL/CL 统计概览 ---
df_cl_sl_summary <- df_qc_a %>%
  summarise(
    Total_Cycles = n(),
    Total_SL_Outlier_Count = sum(Total_SL_Outlier, na.rm = TRUE),
    Total_CL_Outlier_Count = sum(Total_CL_Outlier, na.rm = TRUE),
    CL_AND_SL_Count = sum(CL_AND_SL_Outlier, na.rm = TRUE),
    SL_ONLY_Count = sum(SL_ONLY_Outlier, na.rm = TRUE),
    CL_ONLY_Count = sum(CL_ONLY_Outlier, na.rm = TRUE),
    Total_Final_Outlier = sum(outlier, na.rm = TRUE)
  ) %>%
  mutate(
    Total_SL_Rate = round(Total_SL_Outlier_Count / Total_Cycles * 100, 2),
    Total_CL_Rate = round(Total_CL_Outlier_Count / Total_Cycles * 100, 2),
    CL_AND_SL_Rate = round(CL_AND_SL_Count / Total_Cycles * 100, 2),
    SL_ONLY_Rate = round(SL_ONLY_Count / Total_Cycles * 100, 2),
    CL_ONLY_Rate = round(CL_ONLY_Count / Total_Cycles * 100, 2)
  )

# --- 5.3 提取 Cycle ID 列表 (仅展示前 5 个 ID) ---
get_cycle_ids <- function(flag_col, n=5) {
  ids <- df_qc_a %>%
    filter(!!sym(flag_col)) %>%
    select(Cycle_ID) %>%
    head(n) %>%
    pull()
  return(paste(ids, collapse = ", "))
}

df_cycle_id_analysis <- tribble(
  ~Category, ~Definition, ~Cycle_ID_Example,
  "CL-Only Outliers", "Total_CL_Outlier & !Total_SL_Outlier", get_cycle_ids("CL_ONLY_Outlier"),
  "SL-Only Outliers", "Total_SL_Outlier & !Total_CL_Outlier", get_cycle_ids("SL_ONLY_Outlier"),
  "CL & SL Overlap", "Total_CL_Outlier & Total_SL_Outlier", get_cycle_ids("CL_AND_SL_Outlier")
)


# ----------------------------------------------------------------------
# 6. 输出结果
# ----------------------------------------------------------------------

message("\n--- 阶段 I：QC-A 质控结果 (部分展示) ---")
print(head(df_qc_a %>% select(Cycle_ID, year, Z_Composite_D, CL_D, SL_D, Total_CL_Outlier, Total_SL_Outlier, CL_AND_SL_Outlier, SL_ONLY_Outlier, CL_ONLY_Outlier, outlier), 10))


message("\n--- 统计分析：按年份统计总 Outlier 数量 ---")
print(df_yearly_stats)

message("\n--- 异常 Cycle ID 示例 (前 5 个 Cycle ID) ---")
print(df_cycle_id_analysis)

# 沿用原代码的 SL 独立统计输出
message("\n--- 统计分析：静态规格限 (SL) 独立统计 (总计) ---")
df_sl_stats_total <- df_qc_a %>%
  summarise(
    Total_Cycles = n(),
    SL_D_Outlier_Count = sum(z_conposite_D_SL, na.rm = TRUE),
    SL_O18_Outlier_Count = sum(z_conposite_O18_SL, na.rm = TRUE)
  ) %>%
  mutate(
    SL_D_Rate = round(SL_D_Outlier_Count / Total_Cycles * 100, 2),
    SL_O18_Rate = round(SL_O18_Outlier_Count / Total_Cycles * 100, 2)
  )
print(df_sl_stats_total)

message("\n--- 统计分析：静态规格限 (SL) 独立统计 (按年份) ---")
df_sl_yearly_stats <- df_qc_a %>%
  group_by(year) %>%
  summarise(
    Total_Cycles = n(),
    SL_D_Outlier_Count = sum(z_conposite_D_SL, na.rm = TRUE),
    SL_O18_Outlier_Count = sum(z_conposite_O18_SL, na.rm = TRUE)
  ) %>%
  mutate(
    SL_D_Rate = round(SL_D_Outlier_Count / Total_Cycles * 100, 2),
    SL_O18_Rate = round(SL_O18_Outlier_Count / Total_Cycles * 100, 2)
  )
print(df_sl_yearly_stats)

message("\n--- 统计分析：SL 和 CL 异常的集合分析 (总计) ---")
print(df_cl_sl_summary)