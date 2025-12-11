### Prompt for QC-A
    结合附件“Range-Difference Outlier Detection Script.r”的代码，以及其他的QC相关的附件：
    “Prompt for QC check.txt”
    “/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC_VF-R2-DWD-L1-MD_20251205.md”
    “/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/Range-Difference Outlier Detection Script.r”
    “WVISS.all.csv”
    结合“WVISS.all.csv”和“Prompt for QC check.txt”给我 阶段 I：QC-A质控的R代码，要求代码中包含详细注释说明每一步的功能和目的。
    根据D和18O分别计算SD和Range，计算Z_SD和Z_Range，进而计算Z_Composite。然后根据前五年最稳定数据计算规格限SL和动态控制限CL，最后对每个Cycle进行QC-A质控判定，输出质控结果。分别标记D和18O的质控结果outlier or not，如果一组cycle 的D或18O6任一一种异常，该cycle标记为outlier。
    请确保代码能够完整运行，并生成所需的QC-A质控结果。
    
    我现在出现了一个疑惑，或者反思：如果基于前五年的数据来界定全局异常，那么意味着后五年的大部分cycle将会被识别为outlier。
    我的想法是，后五年的range和sd虽然波动范围大，但是某些cycle之间的HiR是稳定出现的，那么，实际上，这些cycle是可以用于校正的。
    这里面的存在一个理论博弈：理论上，后5年的波动大，意味着一次测量的不确定性大大超过前5年，从概率上讲，cycle之间的HiR稳定出现的可能性应该很小。
    但是，某些时段的确出现了连续稳定的HiR cycles，基于刚才的概率推导，这种偶然性是极低的，但事实确实是出现了连续稳定的HiR cycles,那么，是不是意味着这些cycles连续稳定的HiR是可靠的？
    还是说，这些cycles的HiR稳定出现，只是巧合，是不可靠的？当然，刚才推了一遍这种巧合的概率极低，从统计的角度而言，就可以推断这种稳定不是巧合。因此，反倒说明这种“巧合”是机器的稳定表现，是可以用于WVISS和AWV校正的cycles。
    综上所述，如果只基于稳定期（假设就是前五年）的全局数据来界定异常，那么后五年的大部分cycles将会被识别为异常，从而无法用于校正。
    所以我想讨论的是，测量精度（误差） 和 测量稳定性（HiR稳定出现） 之间的关系。最后，是不是可以通过WVISS建立的HiR校正AWV后D和18O的关系是否能复现 全球降水线GMWL的 8:1关系来反映 测量精度（误差） 和 测量稳定性（HiR稳定出现） 权衡对AWV校正的影响？
    或者说，是不是一定得测量精度（误差）高，才能用于AWV校正；又或许，确实在低精度的情况下，如果确实出现了连续稳定的HiR cycles，也能用于校正，因为这种“巧合的稳定”就不再是巧合，而是一种可靠的稳定？
    这里我要补充一点，因为5个gradient是基于100个观测值所计算的均值，因此每个raw_bias还有一个与之相匹配的精度，即，我所计算的SD，例如一个dO_raw_bias_G5就有一个dO_raw_bias_G5_SD与之对应。
    我自己最近也意识到了，我的R2，DWD，L1，MD等方法，都没有考虑raw_bias 的SD，所以，增加SD的权衡，应该是一个新的视角，就是引进了 “测量精度（误差） 和 测量稳定性（HiR稳定出现）之间的权衡？你怎么看，给我系统分析下。
    其实还有另外一层权衡，严格质控 vs 稳定性HiR下的尽量保留数据。 为什么会有这个权衡，因为后5年的数据确实波动大，HiR不稳定，如果严格质控，几乎就剩不下多少数据可以参与AWV得校正。因为，所以的WVISS质控和HiR建模，都是为了能精准的校正AWV，如果都是outlier，后5年的数据基本就是全废了。
    当然了，最后还得看看在各种质控（宽松、一般、严格的质控）后，到底还能有多少WVISS可以参与AWV的校正。

    cd "/Applications/Working documents/Stable isotopes/江西千烟洲/Atmosphere vapr calibration/QC scripts/"
    pandoc -s "WVISS质控框架的综合分析：精度、稳定性与WLS权衡（完整推理版）.md" -o "WVISS质控框架的综合分析：精度、稳定性与WLS权衡（完整推理版）.docx"


    为了确定变异性增大，是不是可以先计算每一个cycle 的SD和range，然后计算一个median_SD_10year和median_range_10year，画出 range~DOY 和 SD~DOY的10年时间序列，根据与median_SD_10year或median_range_10year的大小关系家颜色；其次，是否需要计算滑动平均趋势，例如，以210个cycle为单位窗口？
    是否还有展示仪器变异性的其他计算指标的建议，以及视图呈现建议？
    先讨论方案，再写代码