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


# ======================================================================
# WVISS QC 参数 (N 和 k) 敏感性分析 R 脚本
# 目标: 
# 1. 确定最小稳定窗口宽度 N (基于 Trend_MAD)
# 2. 确定最大可接受步长 k (基于最大绝对误差 MAE)
# ======================================================================

# 1. 加载所需的库
# 确保已安装这些包: install.packages(c("dplyr", "zoo", "ggplot2", "RcppRoll"))
library(dplyr)
library(zoo)
library(ggplot2)

# 如果 RcppRoll 未安装，将使用 zoo::rollapply，性能会慢一些
tryCatch(library(RcppRoll), error = function(e) message("RcppRoll not found. Using zoo::rollapply (slower) for rolling median."))

# ----------------------------------------------------------------------
# 2. 定义核心统计函数
# ----------------------------------------------------------------------

# 鲁棒 MAD 函数 (按惯例乘以 1.4826 接近正态分布的 SD)
mad_robust <- function(x) {
  # constant = 1.4826, scale to match the standard deviation for normal distribution
  mad(x, constant = 1.4826, na.rm = TRUE)
}

# ----------------------------------------------------------------------
# 3. 数据模拟 (请用您的 WVISS.all.csv 中的实际 Z_Composite 数据替换此块)
# ----------------------------------------------------------------------
set.seed(123)
n_cycles <- 3000 # 假设有 3000 个 Cycle 数据点
# 模拟 Z_Composite 数据：包含一个轻微的线性退化趋势和噪声
df_sim <- data.frame(
  Cycle = 1:n_cycles,
  # Z_Composite 值的模拟 (pmax(0, ...) 确保结果非负)
  Z_Composite_D = pmax(0, rnorm(n_cycles, mean = 0.5 + 0.0005 * 1:n_cycles, sd = 1.5)) 
)
message(paste("已使用", n_cycles, "个 Cycle 模拟数据进行敏感性分析。"))

# 用于分析的变量
DATA_VECTOR <- df_cycle_metrics[["Z_Composite_D"]]

# ======================================================================
# 第一部分：滑动窗口宽度 N 敏感性分析 (寻找 Trend_MAD 拐点)
# ======================================================================

# 4. 定义候选 N 值和计算 Trend_MAD
candidate_N <- seq(1,361,by=10) # 候选窗口大小 (Cycles)
results_N <- data.frame(N = candidate_N, Trend_MAD = NA)

for (i in 1:length(candidate_N)) {
  N_window <- candidate_N[i]
  
  # 计算滚动中位数趋势线 (使用 k=1, 即全分辨率滚动)
  # roll_median/rollapply 的 width 参数是窗口大小 N
  if (exists("roll_median", mode = "function")) {
    rolling_median <- RcppRoll::roll_median(DATA_VECTOR, n = N_window, align = "right", fill = NA)
  } else {
    rolling_median <- zoo::rollapply(DATA_VECTOR, width = N_window, FUN = median, align = "right", fill = NA, na.rm = TRUE)
  }
  
  # 计算 Trend_MAD(N): 趋势线本身的波动性
  results_N$Trend_MAD[i] <- mad_robust(rolling_median)
}

print("--- N 敏感性分析结果 (Trend_MAD) ---")
print(results_N)

head(results_N)
# 5. 可视化 N 敏感性分析结果 (寻找拐点)
plot_N_sensitivity <- ggplot(results_N, aes(x = N, y = Trend_MAD)) +
  geom_point(color = "#0072B2", size = 3) +
  geom_line(color = "#0072B2", size = 1) +
  geom_vline(xintercept = 210, linetype = "dashed", color = "red") + # 假设 210 是选定的 N
  labs(
    title = paste0("窗口宽度 N 敏感性分析 (基于 ", ANALYSIS_COL, ")"),
    subtitle = "Trend_MAD (趋势线波动性) 随 N 变化趋势 - 寻找趋于平稳的拐点",
    x = "滑动窗口宽度 N (Cycles)",
    y = "滚动中位数趋势线的 MAD (Trend_MAD)"
  ) +
  theme_minimal()

# 推荐: 绘制此图表以肉眼识别曲线开始平稳的拐点，该点即为最小稳定 N
print(plot_N_sensitivity)

# ----------------------------------------------------------------------
# 第二部分：滑动步长 k 敏感性分析 (评估滞后误差 MAE)
# ----------------------------------------------------------------------

# 6. 确定基准趋势线 (使用选定的 N 和 k=1)
# 假设根据 N 敏感性分析，我们选择 N_selected
N_selected <- 210 

if (exists("roll_median", mode = "function")) {
  Trend_k1 <- RcppRoll::roll_median(DATA_VECTOR, n = N_selected, align = "right", fill = NA)
} else {
  Trend_k1 <- zoo::rollapply(DATA_VECTOR, width = N_selected, FUN = median, align = "right", fill = NA, na.rm = TRUE)
}

# 7. 定义候选 k 值和计算最大绝对误差 (MAE)
# k=1 是基线 (误差为 0)，k=N_selected 意味着每 N 个 Cycle 才更新一次
candidate_k <- seq(1, N_selected * 2,by=5) 
results_k <- data.frame(k = candidate_k, Max_Absolute_Error = NA)

for (i in 1:length(candidate_k)) {
  k_step <- candidate_k[i]
  
  if (k_step == 1) {
    results_k$Max_Absolute_Error[i] <- 0
    next
  }
  
  # 1. 确定稀疏更新点 (每 k_step 更新一次 CL)
  # 仅在索引 i, i+k, i+2k, ... 处计算滚动中位数
  calculated_indices <- seq(N_selected, length(DATA_VECTOR), by = k_step) 
  
  # 2. 计算稀疏趋势值
  trend_k_sparse <- rep(NA, length(calculated_indices))
  for (j in 1:length(calculated_indices)) {
    idx <- calculated_indices[j]
    start_idx <- idx - N_selected + 1
    window_data <- DATA_VECTOR[start_idx:idx]
    trend_k_sparse[j] <- median(window_data, na.rm = TRUE)
  }
  
  # 3. 提取 k=1 趋势线在相同点的数值
  Trend_k1_sparse_points <- Trend_k1[calculated_indices]
  
  # 4. 计算最大绝对误差 MAE: 量化由于更新滞后造成的最大偏差
  error_abs <- abs(Trend_k1_sparse_points - trend_k_sparse)
  
  results_k$Max_Absolute_Error[i] <- max(error_abs, na.rm = TRUE)
}

print(paste0("--- k 敏感性分析结果 (基于 N=", N_selected, " 的最大绝对误差 MAE) ---"))
print(results_k)


# 假设 Trend_k1 已经计算完成，且已去除 NA
Trend_k1 <- Trend_k1[!is.na(Trend_k1)]
data_length <- length(Trend_k1)

for (i in 1:length(candidate_k)) {
  k_step <- candidate_k[i]
  
  if (k_step == 1) {
    results_k$Max_Absolute_Error[i] <- 0
    next
  }
  
  # 1. 识别最近更新点的索引 i_nearest (关键步骤)
  # i_nearest 是 Cycle i 对应的最近的 k 的倍数 (例如 k=10: Cycle 15 对应 i_nearest=10)
  # seq_len(data_length) 得到的是索引 1, 2, 3, ..., data_length
  
  # i_nearest 向量的计算：
  # (idx - 1) %/% k_step * k_step + 1
  # (例如，idx=15, k_step=10: (15-1)//10 * 10 + 1 = 1*10 + 1 = 11. Cycle 15 的值应追溯到 Cycle 11 的更新)
  
  indices <- seq_len(data_length)
  i_nearest <- floor((indices - 1) / k_step) * k_step + 1
  # 确保 i_nearest 不超过 Trend_k1 的长度 (通常不会，但预防性处理)
  i_nearest[i_nearest > data_length] <- data_length
  
  # 2. 构造稀疏更新趋势线 Trend_sparse_k (与 Trend_k1 长度相同)
  # Trend_k1[i_nearest] 即为：所有 Cycle i 都采用最近更新点 i_nearest 的 Trend_k1 值
  Trend_sparse_k <- Trend_k1[i_nearest]
  
  # 3. 计算最大绝对误差 MAE: 比较两个完整向量
  absolute_error <- abs(Trend_k1 - Trend_sparse_k)
  
  results_k$Max_Absolute_Error[i] <- max(absolute_error, na.rm = TRUE)
}

print(paste0("--- k 敏感性分析结果 (基于 N=", N_selected, " 的最大绝对误差 MAE) ---"))
print(results_k)

quartz("results_k")
plot(Max_Absolute_Error~k,data=results_k)
abline(h=min(results_k$Max_Absolute_Error),col=2)
# 8. 可视化 k 敏感性分析结果 (寻找误差控制阈值)
# 假设可接受的最大误差阈值 (例如 Trend_k1 的 5% MAD)
Error_Threshold <- mad_robust(Trend_k1) * 0.05 

plot_k_sensitivity <- ggplot(results_k, aes(x = k, y = Max_Absolute_Error)) +
  geom_point(color = "#D55E00", size = 3) +
  geom_line(color = "#D55E00", size = 1) +
  geom_hline(yintercept = Error_Threshold, linetype = "dotted", color = "black", size = 1) +
  labs(
    title = paste0("滑动步长 k 敏感性分析 (基于 N=", N_selected, ")"),
    subtitle = paste0("最大绝对误差 (MAE) 随 k 变化趋势。误差阈值 = ", round(Error_Threshold, 4)),
    x = "滑动步长 k (Cycles)",
    y = "最大绝对误差 MAE (相对 k=1 趋势线)"
  ) +
  theme_minimal()
quartz()
# 推荐: 绘制此图表以确定 MAE 曲线不超过误差阈值的最大 k
print(plot_k_sensitivity)

