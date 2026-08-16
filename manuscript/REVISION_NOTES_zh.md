# SJS 重构版说明

## 1. 这次重写完成了什么

本次版本按照 Scandinavian Journal of Statistics 的 Wiley New Journal Design (NJD) 单栏模板重建全文。题目改为 **Debiased Profile Inference for High-Dimensional Mixed-Effects Quantile Regression with Many Small Clusters**，论文身份从“联合惩罚、筛选、算法、CLIME、去偏的模块组合”收束为“many-small-cluster 条件下的高维 profile inference”。

原稿的核心问题已经在公式层面修正：

1. 推断目标改为不依赖平滑带宽的 unsmoothed regularised profile parameter；
2. 样本准则与总体目标统一为等簇权重；
3. 完整 ridge-profile objective 的 score 使用原始固定效应设计 X；
4. profile Hessian 使用严格的 Schur complement；
5. residualised-design 表达补上 A_i^T Lambda A_i；
6. 去偏修正、cluster influence score 和 sandwich variance 全部由同一个 profile objective 推出；
7. CA-IQR-SIS 和 BCD 线性收敛不再作为主要创新；
8. 删除任何依据已知真值校准去偏幅度或标准误的做法；
9. 给出一个明确非空的 bandwidth-sparsity 示例区域。

## 2. 实证数据更换

原 GSE121239 分析被整体替换为 GSE65391。新数据包含 158 名儿科 SLE 患者、924 个 SLE 纵向表达样本，并配有访视、人口学、治疗、疾病活动和连续实验室指标。主响应暂定为 complement C3，主要分析 0.25、0.50 和 0.75 条件分位数，并在数据审计允许时将 0.10 作为敏感性分析；高维预测变量为同次访视的全血转录组表达。

该设计更适合本文：独立簇数量更大；纵向结构清楚；C3 是连续响应，符合局部密度型理论条件；低分位数具有明确的临床解释；临床协变量和治疗信息能够作为不惩罚固定效应纳入模型。

## 3. 当前版本仍未完成的部分

这是一份结构完整、核心公式闭合、且已给出高层条件下主要定理证明的 working manuscript，不是可立即投稿的最终稿。下列内容被有意保留为待完成状态，未虚构数值：

- local profile-Hessian stability、score remainder 与 variance plug-in stability 的 primitive-condition 证明；
- corrected estimator 的参考实现及有限差分单元测试；
- 至少 500/1000 次的完整模拟；
- GSE65391 的缺失数据审计、最终 n/N、模型拟合和结果图表；
- split-sample exploratory gene inference 的聚合规则。

## 4. 推荐的下一阶段验收顺序

第一步先冻结 Proposition 1 的公式并完成数值导数测试。第二步实现 penalised profile estimator、row-wise Dantzig 和 cluster sandwich。第三步用三个 pilot 场景各 200 次重复检查 bias、empirical SD、mean SE 和 coverage。第四步完成 GSE65391 数据审计。第五步再做完整证明审计和全文结果写作。

## 5. 文件说明

- `sjs_profile_quantile_main.tex`：SJS 主文；
- `sjs_profile_quantile_supplement.tex`：补充材料；
- `references.bib`：APA 参考文献库；
- `fetch_wiley_template_files.sh`：获取 Wiley NJDv5 非字体支持文件；
- `lettersp.sty`：Wiley NJD 字距宏兼容文件；
- `compile.sh`：XeLaTeX 编译脚本；
- 编译生成的 PDF 不进入主分支；提交前由源码和冻结结果重新生成。
