## 🔬 QC 参数统计决策分析：滑动窗口 $N$ 与步长 $k$ 的确定
## QC 参数统计决策分析：滑动窗口N与步长k的确定
在科学投稿中，为避免审稿人对滑动窗口参数（$N$ 和 $k$）的主观性质疑，必须进行**参数敏感性分析（Parameter Sensitivity Analysis）**，并用量化统计指标来支撑参数选择。

---

### 1. 窗口宽度 $N$ 的统计确定：趋势线的稳定性

$N$ 的选择目标：找到一个**最小的 $N$**，使得由此计算出的滚动趋势线（$\text{Rolling Median}$）的**波动性达到最小稳定状态**，以确保动态基线（$\text{CL}$）的鲁棒性。

#### 统计指标：趋势线本身的波动性 ($\text{MAD}$)

我们计算滚动中位数趋势线（$\text{Trend}_{\text{full}, 1}$）的**中位数绝对偏差（$\text{MAD}$）**来衡量其波动性。

$$\text{Trend\_MAD}(N) = \text{MAD}(\text{Trend}_{\text{full}, 1} \text{ for window } N)$$

#### 决策标准：拐点法则 (Inflection Point Rule)

1. **方法：** 选取一系列候选 $N$ 值（例如 $\text{70, 140, 210, 300, 365}$ Cycles）进行测试。
2. **绘图：** 绘制 $\text{Trend\_MAD}(N)$ 随 $N$ 变化的曲线。
3. **选择：** 选择曲线开始**趋于平稳**的**拐点**对应的 $N$ 值。这个点代表在确保趋势线稳定性的前提下所能选取的**最小 $N$**，避免了不必要的滞后。

---

### 2. 滑动步长 $k$ 的统计确定：滞后误差控制

$k$ 的选择目标：在不牺牲 $\text{QC}$ 敏感性的前提下，找到一个**最大的 $k$** 来提高计算效率。由于 $\text{k}>1$ 带来的误差是**动态 $\text{CL}$ 更新的滞后性**，我们采用**基于误差控制的方法**，避免引入插值偏差。

#### 统计指标：最大绝对误差 ($\text{MAE}$)

$\text{Error}_{\text{max}, k}$ 量化了由于**延迟更新（步长 $k$）**造成的最大趋势信息滞后。它只比较真实计算点。

$$\text{Error}_{\text{max}, k} = \max_{i} \left( \left| \text{Trend}_{\text{full}, 1}[i] - \text{Trend}_{\text{full}, 1}[i_{\text{nearest}}] \right| \right)$$

其中：
* $\text{Trend}_{\text{full}, 1}[i]$：Cycle $i$ 的全分辨率趋势值（$k=1$ 滚动中位数）。
* $\text{Trend}_{\text{full}, 1}[i_{\text{nearest}}]$：距离 Cycle $i$ 最近的那个**稀疏采样点**（$\text{k}$ 的倍数位置）的**全分辨率**趋势值。

#### 决策标准：误差控制最大化原则

1. **设定阈值：** 定义一个可接受的趋势线滞后误差阈值 $\text{Threshold}_{\text{Error}}$。
   $$\text{Threshold}_{\text{Error}} = \alpha \times \text{MAD}_{\text{Overall}}$$
   （例如，$\alpha = 0.1$。误差阈值应小于仪器自身的长期波动性，确保 $\text{QC}$ 敏感度不受影响。）

2. **选择：** 选取一系列候选 $\text{k}$ 值（例如 $\text{k}=2, 5, 7, 10, 14, 21$）。
3. **最优解：** 找到**满足** $\text{Error}_{\text{max}, k} \le \text{Threshold}_{\text{Error}}$ **的最大 $k$ 值**。

---

### 3. 投稿论证结构总结

| 参数 | 统计依据 | 论证目标 | 常用结论 (示例) |
| :--- | :--- | :--- | :--- |
| **$N$ (窗口宽度)** | **$\text{Trend\_MAD}$ (趋势线波动性)** | 确保 $\text{CL}$ 基线的平滑性和鲁棒性。 | $\text{N}=210$ $\text{Cycles}$ 处于 $\text{Trend\_MAD}$ 曲线达到最小稳定状态的**拐点**。 |
| **$k$ (滑动步长)** | **最大绝对误差 ($\text{MAE}$) 控制** | 在可接受的 $\text{QC}$ 误差范围内实现**最高计算效率**。 | $\text{k}=7$ $\text{Cycles}$ 是满足 $\text{Threshold}_{\text{Error}}$ 约束的**最大步长**。 |

pandoc -s "QC 参数统计决策分析：滑动窗口N与步长k的确定.md" -o "QC 参数统计决策分析：滑动窗口N与步长k的确定.docx"
