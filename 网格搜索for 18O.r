# ----------------------------------------------------------------------
# 步骤 0: 环境准备与参数设置
# ----------------------------------------------------------------------

# 确保 RcppRoll 包已安装和加载，用于高效计算滑动中位数
if (!requireNamespace("RcppRoll", quietly = TRUE)) {
  install.packages("RcppRoll")
}
library(RcppRoll)

# 假设您的数据向量已加载到环境中
# DATA_VECTOR <- ... # 请确保您的实际数据已读入此变量中

DATA_VECTOR <- df_cycle_metrics[["Z_Composite_O18"]]


# 定义综合指标的权重 (侧重异常值检测)
W1_TREND_MAD <- 0.3
W2_MAX_AE <- 0.7

# ----------------------------------------------------------------------
# 步骤 1: 定义参数空间 (N 和 k)
# ----------------------------------------------------------------------

# N 窗口长度：从 1 到 401，步长 10
N_list <- seq(7, 365, by = 5)
# k 步长：从 1 到 N_max (401)，步长 5
k_list_base <- seq(2, 365, by = 5)

# 初始化结果存储 (包含所有 N 和 k 组合)
Results_S <- data.frame(
  N = integer(),
  k = integer(),
  Trend_MAD = numeric(),
  Max_AE = numeric(),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------
# 步骤 2: 二维网格搜索 (计算 Trend_MAD 和 Max_AE)
# ----------------------------------------------------------------------

L <- length(DATA_VECTOR)

for (N in N_list) {
  # 计算 N 窗口长度下的滑动中位数 (k=1 全分辨率基线)
  # 使用 align = "right" 实现因果式滑动窗口 (Trailing Moving Median)
  rolling_median_k1 <- roll_median(DATA_VECTOR, n = N, align = "right", fill = NA)
  
  # 遍历 k 列表
  for (k in k_list_base) {
    # 约束条件: k 不能超过 N
    if (k > N) {
      next
    }
    
    # ----------------------------------------------------------
    # 核心计算逻辑 (基于 k > 1 的稀疏更新点)
    # ----------------------------------------------------------
    
    # 确定稀疏更新点 (仅在 i = N, N+k, N+2k, ... 处计算残差)
    calculated_indices <- seq(N, L, by = k)
    
    # 提取实际观测值 Y_i (在更新点处)
    Y_i_sparse <- DATA_VECTOR[calculated_indices]
    
    # 提取趋势估计值 hat{Y}_i (在更新点处)
    # hat{Y}_i 在 k>1 时仍使用 N 个点计算，但只在稀疏点上采样
    hat_Y_i_sparse <- rolling_median_k1[calculated_indices]
    
    # 确保没有 NA 值 (rolling_median_k1 在 i < N 时是 NA)
    valid_indices <- !is.na(hat_Y_i_sparse)
    Y_i_sparse <- Y_i_sparse[valid_indices]
    hat_Y_i_sparse <- hat_Y_i_sparse[valid_indices]
    
    # 如果没有有效点，则跳过
    if (length(Y_i_sparse) == 0) {
      next
    }
    
    # 计算稀疏残差序列
    residuals_sparse <- abs(Y_i_sparse - hat_Y_i_sparse)
    
    # 3.2. 计算指标
    Current_Trend_MAD <- median(residuals_sparse, na.rm = TRUE)
    Current_Max_AE  <- max(residuals_sparse, na.rm = TRUE)
    
    # 3.3. 存储原始结果
    Results_S <- rbind(Results_S, data.frame(
      N = N,
      k = k,
      Trend_MAD = Current_Trend_MAD,
      Max_AE = Current_Max_AE
    ))
  }
}

# ----------------------------------------------------------------------
# 步骤 3: 归一化和计算综合指标 S(N, k)
# ----------------------------------------------------------------------

# 3.1. 归一化 (找到全局 Min/Max)
X_min_TMAD <- min(Results_S$Trend_MAD)
X_max_TMAD <- max(Results_S$Trend_MAD)
X_min_MAE <- min(Results_S$Max_AE)
X_max_MAE <- max(Results_S$Max_AE)

# 3.2. 计算归一化值
Results_S$Trend_MAD_norm <- (Results_S$Trend_MAD - X_min_TMAD) / (X_max_TMAD - X_min_TMAD)
Results_S$Max_AE_norm <- (Results_S$Max_AE - X_min_MAE) / (X_max_MAE - X_min_MAE)

# 3.3. 计算综合指标 S(N, k)
Results_S$S_score <- (W1_TREND_MAD * Results_S$Trend_MAD_norm) + (W2_MAX_AE * Results_S$Max_AE_norm)

# ----------------------------------------------------------------------
# 步骤 4: 确定最优解
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# 步骤 4: 确定最优解与鲁棒解集
# ----------------------------------------------------------------------

# 找到 S 分数的最小值
S_min <- min(Results_S$S_score)
Optimal_Result <- Results_S[which.min(Results_S$S_score), ]

# 计算 S 分数的标准差 (衡量解的集中度/波动性)
SD_S <- sd(Results_S$S_score)

# 定义鲁棒性阈值: S_threshold = S_min + C * SD_S
# C 是一个经验常数，通常取 0.5 到 1.5。这里我们取 C=1 (S_min + 1*SD)
C_ROBUSTNESS <- 1 
S_threshold <- S_min + C_ROBUSTNESS * SD_S

# 筛选出鲁棒解集 (所有 S_score 低于阈值的组合)
Robust_Solutions <- Results_S[Results_S$S_score <= S_threshold, ]

# ----------------------------------------------------------------------
# 步骤 5: 结果分析与输出
# ----------------------------------------------------------------------

cat("--- 滑动窗口双参数 (N, k) 综合优化结果 (包含鲁棒性筛选) ---\n")
cat(sprintf("权重设置: Trend_MAD (w1=%.1f), Max_AE (w2=%.1f)\n", W1_TREND_MAD, W2_MAX_AE))
cat(sprintf("综合 S 分数标准差 (SD_S): %.4f\n", SD_S))
cat(sprintf("最优 S 分数 (S_min): %.4f\n", S_min))
cat(sprintf("鲁棒性阈值 (S_min + 1*SD_S): %.4f\n\n", S_threshold))

cat("1. 全局最优解 (S_min):\n")
print(Optimal_Result)

cat("\n2. 鲁棒解集 (S_score <= 阈值) 摘要:\n")

# 计算鲁棒解集的关键统计量
N_robust_min <- min(Robust_Solutions$N)
N_robust_max <- max(Robust_Solutions$N)
k_robust_min <- min(Robust_Solutions$k)
k_robust_max <- max(Robust_Solutions$k)

cat(sprintf("  - N 鲁棒范围: [%d, %d]\n", N_robust_min, N_robust_max))
cat(sprintf("  - k 鲁棒范围: [%d, %d]\n", k_robust_min, k_robust_max))
cat(sprintf("  - 解集大小: %d 个组合 (占总 %d 个组合的 %.2f%%)\n", 
            nrow(Robust_Solutions), nrow(Results_S), 
            nrow(Robust_Solutions) / nrow(Results_S) * 100))

# 推荐的最佳鲁棒解 (例如: 鲁棒集中 k 最小、N 适中的点)
# 推荐策略: 在鲁棒解集中，选择 k 最小，同时 N 适中（靠近平均或中位数）的解
N_median_robust <- median(Robust_Solutions$N)
Recommended_Solution <- Robust_Solutions[
  order(Robust_Solutions$k, abs(Robust_Solutions$N - N_median_robust)), 
][1, ]

cat("\n3. 推荐的鲁棒最优解:\n")
print(Recommended_Solution)


#最优解热力图
# ----------------------------------------------------------------------
# 步骤 1: 数据准备 - 将 S_score 转换为矩阵 (Data Preparation)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# 步骤 1: 数据准备 - 将 S_score 转换为矩阵 (Data Preparation)
# ----------------------------------------------------------------------

# 确保 N_unique 和 k_unique 已经定义
N_unique <- unique(Results_S$N)
k_unique <- unique(Results_S$k)

# 重建 S_score 矩阵
S_matrix <- matrix(NA, nrow = length(N_unique), ncol = length(k_unique))

for (i in 1:nrow(Results_S)) {
  N_idx <- which(N_unique == Results_S$N[i])
  k_idx <- which(k_unique == Results_S$k[i])
  S_matrix[N_idx, k_idx] <- Results_S$S_score[i]
}

# ----------------------------------------------------------------------
# 步骤 2: 自定义颜色梯度和分界点
# ----------------------------------------------------------------------

# 2.1 定义 S_score 范围的最大值
S_max <- max(Results_S$S_score)

# 2.2 定义非鲁棒区域的梯度数量 (取 8 个)
N_GRADIENTS_REST <- 8

# 2.3 定义分界点 (Breaks)
# 第一个分界点: S_min
# 第二个分界点: S_threshold (S_min + SD_S)
breaks <- c(
  S_min, 
  S_threshold, 
  # 剩下的 N_GRADIENTS_REST 个分界点
  seq(S_threshold + (S_max - S_threshold) / N_GRADIENTS_REST, S_max, length.out = N_GRADIENTS_REST)
)
# 确保分界点唯一且排序
breaks <- sort(unique(breaks))

# 2.4 定义颜色 (Colors, 数量 = length(breaks) - 1)
# 颜色 1: 鲁棒解集 (S <= S_threshold) - 固定的深蓝色 (最优)
color_robust <- "#1A5276" 
# 颜色 2 到 9: 非鲁棒解集 - 从浅蓝过渡到红 (性能逐渐变差)
color_rest_palette <- colorRampPalette(c("#A9CCE3", "yellow", "red"))(length(breaks) - 2)
custom_colors <- c(color_robust, color_rest_palette)

# ----------------------------------------------------------------------
# 步骤 3: 使用 layout() 绘制热力图和颜色图例
# ----------------------------------------------------------------------

quartz("heatmap of solutions for 18O S-score")
# 设置布局: 5单位宽度给图，1单位宽度给颜色图例
layout(matrix(c(1, 2), ncol = 2), widths = c(5, 1))
# 3.1 绘制主图 (N-k 平面上的 S_score)
par(mar = c(5, 5, 4, 1) + 0.1) # 调整右侧边距

# !!! 修复错误: 移除 useRaster = TRUE !!!
image(
  x = N_unique, 
  y = k_unique, 
  z = S_matrix, 
  breaks = breaks,        # 使用自定义分界点
  col = custom_colors,    # 使用自定义颜色
  xlab = "Window Length N (Cycles)", 
  ylab = "Step Size k (Cycles)",
  main = "Contour Map of Composite Score S(N, k)"
)

# 3.2 突出显示鲁棒解集和全局最优解
Robust_N <- Results_S$N[Results_S$S_score <= S_threshold]
Robust_k <- Results_S$k[Results_S$S_score <= S_threshold]

# 鲁棒解集 (半透明白色圆点)
points(
  x = Robust_N, 
  y = Robust_k, 
  pch = 16, 
  col = rgb(1, 1, 1, 0.5), 
  cex = 0.8
)

# 全局最优解 (红色叉号)
points(
  x = Optimal_Result$N, 
  y = Optimal_Result$k, 
  pch = 4, 
  col = "red", 
  cex = 2, 
  lwd = 3
)

# ----------------------------------------------------------------------
# 步骤 3.3 绘制颜色图例 (Color Bar) - 使用 rect() 绘制
# ----------------------------------------------------------------------

par(mar = c(5, 0, 4, 3) + 0.1) # 调整边距 (左侧边距为 0)

# 创建一个空的绘图区域来承载图例
plot.new()
# 明确设置用户坐标系，这是正确的步骤
plot.window(xlim = c(0, 1), ylim = range(breaks))

# 使用 rect() 绘制颜色条
# breaks 定义了颜色段的边界
for (i in 1:(length(breaks) - 1)) {
  rect(
    xleft = 0,                # 左边界
    ybottom = breaks[i],      # 当前段的下边界
    xright = 1,               # 右边界
    ytop = breaks[i+1],       # 当前段的上边界
    col = custom_colors[i],   # 对应颜色
    border = NA               # 无边框
  )
}

# 添加轴标签
axis(4, at = breaks, labels = format(breaks, digits = 3), las = 1) # las=1 使标签水平显示
mtext("S Score", side = 4, line = 4, cex = 1)
box()

# 添加鲁棒阈值标记 (虚线)
abline(h = S_threshold, col = "black", lwd = 2, lty = 2)
mtext("Robust Threshold", side = 4, line = 3, at = S_threshold, col = "black", cex = 0.8)

# 恢复默认布局
layout(matrix(1), widths = 1, heights = 1)
