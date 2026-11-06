# Loop 入组两周时的夜间血糖预测

探索侧登记 404 人，达到监测条件的患者有 252 人，每人提供三个时段，共 756 个结局。216 名患者留作验证。

到第 14 日结束时，以已经取得的监测数据预测随后三个 21 日时段的夜间高值。每夜计算最高 10% 血糖读数的均值，每个时段对合格夜间等权平均。评价以患者等权的第 75 百分位预测损失为准，数值越小越好。

CQ_SMOOTH 的损失为 6.7914。临床分位数回归为 6.8248，经分位数校准的 elastic net 为 7.3754。相对改善分别为 +0.49% 和 +7.92%。

每折训练结束后，288 个血糖曲线变量中系数绝对值超过 10⁻⁷ 的个数依次为 1, 0, 0, 0, 0。Python 与 R 分别重建了预测误差，DuckDB 从五分钟时格重新计算了夜间结局和 288 个输入变量。

必要比较未达到登记的 2% 条件。保留患者的随访结局继续封存，本次计算未形成经独立验证的临床发现。

为核查完整曲线能否提供额外的预测信息，在相同患者和相同训练组中另行拟合了只使用临床变量的校准 elastic net。其损失为 7.0562，加入 288 个曲线变量后为 7.3754，相对改善为 -4.52%。

Lum 等，Diabetes Technology & Therapeutics，2021，DOI 10.1089/dia.2020.0535。

The source of the data is the Loop Study (sponsored by the Jaeb Center for Health Research and funded by the Helmsley Charitable Trust), but the analyses, content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the study sponsor.
