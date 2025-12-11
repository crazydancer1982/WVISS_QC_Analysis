# ----------------------------------------------------------------------
# 阶段 I：QC-A (VF - 仪器波动与绝对合格性质控) R 代码
# 目标: 结合 WVISS.all.csv 真实数据，计算 Z_Composite，确定 SL/CL 阈值，
#       并对每个 Cycle 进行 QC-A 判定。
# ----------------------------------------------------------------------

# 加载所需的库
# 确保已安装：dplyr, tidyr, zoo
library(dplyr)
library(tidyr)
library(zoo)

options("digits" = 9) 

# ----------------------------------------------------------------------
# 1. 实际数据加载与预处理 (替换模拟数据块)
# ----------------------------------------------------------------------

# **请手动将此路径修改为 WVISS.all.csv 所在的实际目录**
wd.synthesis <- '/Applications/Working documents/Stable isotopes/江西千烟洲/synthesized datasets_2011-20'
setwd(wd.synthesis)

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

# --- 1.1 核心偏差计算 (沿用 Shiny App 逻辑) ---
# 定义 SMOW 比值 (来自 Shiny App)
RD_smow <- 0.00015576
RO_smow <- 0.0020052

WVISS.all <- WVISS.all %>%
  # 计算观测值 dO_obs 和 dD_obs
  mutate(
    dO_obs = ((O18ppm_mean / O16ppm_mean) / RO_smow - 1) * 1000,
    dD_obs = ((HODppm_mean / (O16ppm_mean * 2)) / RD_smow - 1) * 1000,
    # 计算原始偏差 (Observed - True)
    dO_raw_bias = dO_obs - o_wahaha,
    dD_raw_bias = dD_obs - d_wahaha
  ) %>%
  # 确保 group.id 是唯一 Cycle 标识符，Gradient 是 Cycle 内的 1-5 标气点
  select(group.id, Gradient, dD_raw_bias, dO_raw_bias, year) 

# --- 1.2 数据重塑：长格式转宽格式 ---
# 目标：将 Cycle 内的 5 个 Gradient 偏差值 (dD_raw_bias, dO_raw_bias) 转换成单行数据，
#       以便于计算 Cycle 内部的 SD 和 Range。
df_wide <- WVISS.all %>%
  pivot_wider(
    id_cols = c(group.id, year),
    names_from = Gradient,
    values_from = c(dD_raw_bias, dO_raw_bias),
    names_prefix = "G"
  ) %>%
  # 确保 Cycle 顺序正确 (按时间或 group.id 排序)
  arrange(group.id)

# ----------------------------------------------------------------------
# 2. 准备 QC-A 核心指标数据：SD 和 Range
# ----------------------------------------------------------------------

# 定义同位素类型和相应的 G 列名
isotopes <- c("D", "O18")
G_cols <- list(
  D = paste0("dD_raw_bias_G", 1:5),
  O18 = paste0("dO_raw_bias_G", 1:5)
)

# 初始化结果数据框 (以宽格式 df_wide 为基础)
df_qc_a <- df_wide %>% rename(Cycle_ID = group.id)

for (iso in isotopes) {
  # (a) 计算 Cycle 内的标准差 (SD)
  # 测量 Cycle 内部数据的平均波动量。
  df_qc_a[[paste0("SD_", iso)]] <- apply(df_qc_a[, G_cols[[iso]]], 1, sd, na.rm = TRUE)
  
  # (b) 计算 Cycle 内的 Range (极差)
  # 测量 Cycle 内部数据的最大绝对波动量 (|G5 - G1|)，沿用 Range-Difference 逻辑。
  df_qc_a[[paste0("Range_", iso)]] <- abs(df_qc_a[[paste0("dD_raw_bias_G", 5)]] - df_qc_a[[paste0("dD_raw_bias_G", 1)]])
  if (iso == "O18") {
    df_qc_a[[paste0("Range_", iso)]] <- abs(df_qc_a[[paste0("dO_raw_bias_G", 5)]] - df_qc_a[[paste0("dO_raw_bias_G", 1)]])
  }
}
str(df_qc_a)
# ----------------------------------------------------------------------
# 3. 计算 Z_SD, Z_Range 和 Z_Composite (波动综合指数)
# ----------------------------------------------------------------------

# **鲁棒性稳定期定义：** 假设前五年为稳定期 (请根据实际情况调整年份)
STABLE_YEAR_CUTOFF <- min(df_qc_a$year) + 4 # 假设前 5 年为稳定期
df_stable <- subset(df_qc_a, year <= STABLE_YEAR_CUTOFF)

# 定义权重：w2 > w1，优先惩罚 Range (极端波动)
w1 <- 0.4 # SD 权重
w2 <- 0.6 # Range 权重

for (iso in isotopes) {
  SD_col <- paste0("SD_", iso)
  Range_col <- paste0("Range_", iso)
  
  # 鲁棒性统计量：计算稳定期的 Median 和 MAD (作为 SD 估算)
  # 使用 MAD (Median Absolute Deviation) 作为比 SD 更稳健的变异性估算。
  mu_sd_stable <- median(df_stable[[SD_col]], na.rm = TRUE)
  sigma_sd_stable <- mad(df_stable[[SD_col]], na.rm = TRUE)
  
  mu_range_stable <- median(df_stable[[Range_col]], na.rm = TRUE)
  sigma_range_stable <- mad(df_stable[[Range_col]], na.rm = TRUE)
  
  # (a) 计算标准化 Z 分数 (基于稳定期鲁棒性参数)
  df_qc_a[[paste0("Z_SD_", iso)]] <- (df_qc_a[[SD_col]] - mu_sd_stable) / sigma_sd_stable
  df_qc_a[[paste0("Z_Range_", iso)]] <- (df_qc_a[[Range_col]] - mu_range_stable) / sigma_range_stable
  
  # (b) 计算 Z_Composite (综合波动指数)
  df_qc_a[[paste0("Z_Composite_", iso)]] <- 
    w1 * df_qc_a[[paste0("Z_SD_", iso)]] + w2 * df_qc_a[[paste0("Z_Range_", iso)]]
}

# ----------------------------------------------------------------------
# 4. 阈值设定：规格限 (SL)
# ----------------------------------------------------------------------

for (iso in isotopes) {
  Z_Comp_col <- paste0("Z_Composite_", iso)
  
  # SL (规格限) X_Z： 基于稳定期 Z_Composite 的 99.73% 分位数 (鲁棒性判定)
  SL_XZ <- quantile(df_stable[[Z_Comp_col]], probs = 0.9973, na.rm = TRUE)
  
  # 存储 SL 值
  df_qc_a[[paste0("SL_XZ_", iso)]] <- SL_XZ
}

# ----------------------------------------------------------------------
# 5. 阈值设定：动态控制限 (CL) - 基于自适应 k 值策略
# ----------------------------------------------------------------------

# 5.1 参数设定
m <- 50 # 滑动窗口大小 m (经验值 30-60 个 Cycle)
k_low_sigma <- 2.5 # 低波动期 k 值 (低误报率)
k_high_sigma <- 1.5 # 高波动期 k 值 (高敏感度)
sigma_threshold_factor <- 1.5 # 切换阈值：当前 MAD > 1.5 * 稳定期 MAD

for (iso in isotopes) {
  Z_Comp_col <- paste0("Z_Composite_", iso)
  
  # 稳定期 Z_Composite 的 MAD
  mad_stable_z <- mad(df_stable[[Z_Comp_col]], na.rm = TRUE)
  
  # 初始化动态 CL 和 k 值
  df_qc_a[[paste0("CL_", iso)]] <- NA
  
  # 循环进行滑动窗口计算和 k 值切换
  for (i in m:nrow(df_qc_a)) {
    window_data <- df_qc_a[(i - m + 1):i, Z_Comp_col]
    
    # 滑动窗口的鲁棒性变异性估算 (使用 MAD)
    mad_current <- mad(window_data, na.rm = TRUE)
    
    # 动态 k 值决策 (自适应敏感度)
    if (mad_current > sigma_threshold_factor * mad_stable_z) {
      k_val <- k_high_sigma # 过程波动大，切换到高敏感度
    } else {
      k_val <- k_low_sigma # 过程稳定，切换到低误报率
    }
    
    # CL 阈值计算 (中心线为 0，因为 Z 分数是标准化后的)
    CL_current <- 0 + k_val * mad_current 
    
    # 赋值给当前周期
    df_qc_a[[paste0("CL_", iso)]][i] <- CL_current
  }
}

# ----------------------------------------------------------------------
# 6. QC-A 质控判定
# ----------------------------------------------------------------------

# 初始化判定结果列
df_qc_a$QC_A_D_Outlier <- FALSE
df_qc_a$QC_A_O18_Outlier <- FALSE
df_qc_a$QC_A_Cycle_Outlier <- FALSE

for (iso in isotopes) {
  Z_Comp_col <- paste0("Z_Composite_", iso)
  SL_col <- paste0("SL_XZ_", iso)
  CL_col <- paste0("CL_", iso)
  
  # 判定逻辑：Z_Composite 必须在动态 CL 内，且不能超过绝对 SL
  df_qc_a[[paste0("QC_A_", iso, "_Outlier")]] <- 
    # 判定 1: 超出绝对规格限 (SL)
    (df_qc_a[[Z_Comp_col]] > df_qc_a[[SL_col]]) | 
    # 判定 2: 超出动态控制限 (CL) - 只对滑动窗口后的数据进行判定
    (df_qc_a[[Z_Comp_col]] > df_qc_a[[CL_col]])
}

# 最终 Cycle 判定：如果一组 Cycle 的 D 或 18O 任一一种异常，该 Cycle 标记为 Outlier
df_qc_a$QC_A_Cycle_Outlier <- df_qc_a$QC_A_D_Outlier | df_qc_a$QC_A_O18_Outlier

# ----------------------------------------------------------------------
# 7. 输出质控结果
# ----------------------------------------------------------------------

# 筛选并输出包含所有 QC 结果的表格的前 10 行和后 10 行
print("--- 阶段 I：QC-A 质控结果 (部分输出) ---")
print(df_qc_a %>%
        select(Cycle_ID, year, starts_with("Z_Composite"), 
               starts_with("SL_XZ"), starts_with("CL_"), 
               starts_with("QC_A_")) %>%
        head(10))

print(df_qc_a %>%
        select(Cycle_ID, year, starts_with("Z_Composite"), 
               starts_with("SL_XZ"), starts_with("CL_"), 
               starts_with("QC_A_")) %>%
        tail(10))

# 总结异常 Cycle 数量
total_outliers <- sum(df_qc_a$QC_A_Cycle_Outlier, na.rm = TRUE)
print(paste("总 Cycles 数:", nrow(df_qc_a)))
print(paste("QC-A 判定异常 Cycles 数:", total_outliers))

# ----------------------------------------------------------------------