# IOBP2 随机分组两周时的夜间血糖预测

探索侧登记 212 人，达到监测条件的患者有 160 人，每人提供三个时段，共 480 个结局。114 名患者留作验证。

到第 14 日结束时，以已经取得的监测数据预测随后三个 21 日时段的夜间高值。每夜计算最高 10% 血糖读数的均值，每个时段对合格夜间等权平均。评价以患者等权的第 75 百分位预测损失为准，数值越小越好。

CQ_SMOOTH 的损失为 7.4488。临床分位数回归为 7.6262，经分位数校准的 elastic net 为 7.1041。相对改善分别为 +2.33% 和 -4.85%。

每折训练结束后，288 个血糖曲线变量中系数绝对值超过 10⁻⁷ 的个数依次为 0, 0, 2, 0, 0。Python 与 R 分别重建了预测误差，DuckDB 从五分钟时格重新计算了夜间结局和 288 个输入变量。

必要比较未达到登记的 2% 条件。保留患者的随访结局继续封存，本次计算未形成经独立验证的临床发现。

为核查完整曲线能否提供额外的预测信息，在相同患者和相同训练组中另行拟合了只使用临床变量的校准 elastic net。其损失为 6.9969，加入 288 个曲线变量后为 7.1041，相对改善为 -1.53%。

Bionic Pancreas Research Group，New England Journal of Medicine，2022，DOI 10.1056/NEJMoa2205225。

The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the analyses, content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the Bionic Pancreas Research Group or Beta Bionics.
