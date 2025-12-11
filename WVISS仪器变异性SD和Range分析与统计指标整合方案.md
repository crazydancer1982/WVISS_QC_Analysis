# 🔬 WVISS 仪器变异性（SD/Range）分析与统计指标整合方案

## 概览：QC 框架中 SL 与 CL 的哲学定位

本方案整合鲁棒性指标计算与可视化，旨在评估 $\text{WVISS}$ 仪器在 10 年间的性能退化模式，并明确 $\text{QC}$ 框架中 **规格限 ($\text{SL}$) ** 和 **控制限 ($\text{CL}$) ** 的哲学定位。

---

## 🛠️ 统计指标计算公式与定义

### 1. 核心鲁棒统计量

| 名称 | 符号 | 计算公式 | 作用 |
| :--- | :--- | :--- | :--- |
| **中位数** | $\tilde{X}$ | $\text{Median}(X_i)$ | 鲁棒的集中趋势，不受极端异常值影响。 |
| **中位数绝对偏差** | $\text{MAD}$ | $\text{Median}(|X_i - \tilde{X}|) \times 1.4826$ | 鲁棒的离散程度衡量，等效于 $\sigma$。 |

### 2. 鲁棒控制限 ($\text{CL}$)

我们使用 **$\text{Median} + 3 \times \text{MAD}$** 作为鲁棒的 $3\sigma$ 控制限（$\text{CL}$），用于**动态过程监控**：

$$\text{Threshold}_{\text{CL}} = \text{Median} + 3 \times \text{MAD}$$

### 3. 综合 $\text{Z}$ 分数 ($\text{Z}_{\text{Composite}}$)

$\text{Z}_{\text{Composite}}$ 是 $\text{QC-A}$ 的**最终决策指标**，用于统一量化 $\text{SD}$ 和 $\text{Range}$ 两个指标的波动程度。

1.  **标准化 $\text{Z}$ 分数（基于 10 年鲁棒统计量）：**

    $$\text{Z}_{\text{SD}, i} = \frac{\text{SD}_{i} - \text{Median}_{\text{SD}}}{\text{MAD}_{\text{SD}}}$$
    $$\text{Z}_{\text{Range}, i} = \frac{\text{Range}_{i} - \text{Median}_{\text{Range}}}{\text{MAD}_{\text{Range}}}$$

2.  **综合指标（取最大绝对值）：**

    $$\text{Z}_{\text{Composite}, i} = \max(|\text{Z}_{\text{SD}, i}|, |\text{Z}_{\text{Range}, i}|)$$
    * **$\text{QC}$ 判定：** 若 $\text{Z}_{\text{Composite}, i} > 3$（或 $3.5$），则判定该 $\text{Cycle}$ 质量失控。

---

## 📈 变异性可视化与趋势分析方案

本方案采用**时间序列点图**、**滑动中位数趋势**和**超限比例分析**的三维视图。

### 1. 核心指标与鲁棒性基线计算

| 步骤 | 计算指标 | 目的 |
| :--- | :--- | :--- |
| **单 $\text{Cycle}$ 指标** | $\text{SD}_{\text{cycle}, i}$ 和 $\text{Range}_{\text{cycle}, i}$ | 衡量第 $i$ 个 $\text{Cycle}$ 的内部波动。 |
| **十年鲁棒基线** | $\text{median\_SD}_{\text{10year}}$ 和 $\text{median\_Range}_{\text{10year}}$ | 计算十年数据的**中位数**作为**鲁棒性**的性能基线。 |

### 2. 时间序列可视化与颜色编码

| 维度 | 描述 | 目的 |
| :--- | :--- | :--- |
| **X 轴** | **日期时间** | 展现时间序列趋势。 |
| **Y 轴** | $\text{SD}_{\text{cycle}, i}$ 或 $\text{Range}_{\text{cycle}, i}$ | 展现变异性大小。 |
| **颜色编码** | **二元编码：** 正常 ($\le \text{Threshold}_{\text{CL}}$) / 警示 ($> \text{Threshold}_{\text{CL}}$) | **即时定位**哪些 $\text{Cycle}$ 的变异性**高于鲁棒控制限**。 |
| **趋势线** | **滑动中位数 ($\text{Rolling Median}$) ** | 鲁棒地追踪仪器**当前**的变异性趋势，用于动态过程监控 ($\text{CL}$)。 |

### 3. 仪器变异性补充指标与视图建议

| 指标/视图 | 建议计算方法 | 分析价值 |
| :--- | :--- | :--- |
| **滑动 $\text{MAD}$ of $\text{SD}$** | 以 $N$ $\text{Cycles}$ 为窗口，计算窗口内 $\text{SD}$ 值的 $\text{MAD}$。 | 衡量仪器**波动性本身的波动性**（$\text{Volatility of Volatility}$），用于判断仪器状态是否趋于不稳定。 |
| **$\text{Z}_{\text{Composite}}$ 质控图** | 绘制 $\text{Z}_{\text{Composite}}$ 时间序列，并标出 $\text{SL}$ 和 $\text{CL}$。 | 最标准的统计控制图，直接用于**硬性排除**或**过程监控**。 |
| **月度/年度超限比例** | 统计每个月/年 $\text{Z}_{\text{Composite}}$ 超出 $\text{Threshold}_{\text{CL}}$ 的 $\text{Cycle}$ 占总 $\text{Cycle}$ 的百分比。 | 宏观评估仪器性能衰退的**持续时间和严重性**。 |

---

## 📋 宏观 $\text{QC}$ 诊断与可视化策略

### 4. 月度/年度精度状态判定 (宏观 QC)

**核心方法：基于 $\text{CL3MAD}$ 的超限百分比**

在 $\text{R}$ 代码中，我们采用以下鲁棒阈值作为判断 $\text{Cycle}$ 质量是否合格的依据：

$$\text{Threshold}_{\text{CL}} = \text{Median} + 3 \times \text{MAD}$$

**判定流程：**

* **Cycle 级判定：** 如果任意同位素（$\text{dD}$ 或 $\text{d}^{18}\text{O}$）的 $\text{SD}$ 或 $\text{Range}$ 超出其各自的 $\text{Threshold}_{\text{CL}}$，则该 $\text{Cycle}$ 被标记为**不合格** ($\text{Outlier}$)。
* **时间窗汇总：** 将 $\text{Outlier}$ $\text{Cycle}$ 的数量按月份或年份分组。
* **百分比计算：**
    $$\text{月度超限比例} (\%\text{Exceed}) = \frac{\text{该月不合格 Cycle 数量}}{\text{该月总 Cycle 数量}} \times 100$$

**诊断价值与管理阈值：**

| 状态 | 比例范围 | 诊断描述 | 建议管理操作 |
| :--- | :--- | :--- | :--- |
| **绿色/稳定** | $0\% \sim 5\%$ | 过程处于统计控制状态，偶然的 $\text{Outlier}$ 是预期内的。 | 持续监控。 |
| **黄色/警示** | $5\% \sim 10\%$ | 过程开始出现**波动加剧**的迹象，性能基础系统性下降。 | 关注并安排维护。 |
| **红色/失控** | $> 10\%$ | 仪器已进入**高风险运行阶段**，数据可用性极低。 | 必须进行维护或校准。 |

### 5. 视觉呈现策略分析

**最佳方案：** **主要展示宏观的月度超限比例（点 $2$），但通过图表设计和注释，间接体现微观 $\text{QC}$ 判定的基础。**

| 图表类型 | 焦点 | 目的 | 视觉元素 |
| :--- | :--- | :--- | :--- |
| **核心图：月度超限比例时间序列** | **宏观风险诊断** | **管理决策**：直接显示仪器精度**系统性失控**的严重性。 | 使用**水平虚线**划分 $5\%$ 和 $10\%$ 管理阈值；使用**颜色编码** (绿/黄/红) 标记数据点状态。 |
| **辅助图：$\text{Z}_{\text{Composite}}$ 质控图** | **微观 $\text{QC}$ 判定** | **验证与追溯**：验证 $\text{Cycle}$ 级判定是否合理，并追溯高超限月份的异常源。 | 绘制 $\text{Z}_{\text{Composite}}$ 散点图，并标出 $\mathbf{Z=3}$ 控制限。 |

---

## 🧭 $\text{SL}$ 与 $\text{CL}$ 的哲学解析与应用

### $\text{SL}$ 与 $\text{CL}$ 的全称和统计哲学

| 缩写 | 英文全称 | 中文译名 | 统计哲学（核心区别） |
| :--- | :--- | :--- | :--- |
| **$\text{SL}$** | **Specification Limit** | **规格限 / 容差限** | 衡量**过程能力**：是否达到了**外部设定的绝对质量标准**。与过程当前的波动无关。 |
| **$\text{CL}$** | **Control Limit** | **控制限** | 衡量**过程稳定性**：是否处于**统计控制状态**。由过程自身的统计波动决定。 |

### $\text{QC}$ 目标与 $\text{SL/CL}$ 的完美挂钩

您的 $\text{QC}$ 策略是一种**先进的 $\text{SPC}$ 实践**，结合了固定 $\text{SL}$ 和动态 $\text{CL}$ 的双重阈值体系：

| 您的 $\text{QC}$ 概念 | 对应的 $\text{SL/CL}$ 哲学 | 作用与匹配性 |
| :--- | :--- | :--- |
| **绝对质控（10年，固定）** | **$\text{SL}$：绝对可接受的规格限** | $\text{SL}$ 基于**仪器历史上的最佳性能**设定**质量底线**。用于**硬性排除**任何达不到基本精度要求的 $\text{Cycle}$。 |
| **滑窗质控（短期，动态）** | **$\text{CL}$：过程自身的控制限** | 滑动 $\text{CL}$ 持续追踪仪器**当前**的变异性。用于**警示**：如果一个 $\text{Cycle}$ 异于其邻近周期的平均波动，则表明仪器状态**发生突变**。 |

pandoc -s "WVISS仪器变异性SD和Range分析与统计指标整合方案.md" -o "WVISS仪器变异性SD和Range分析与统计指标整合方案.docx"

