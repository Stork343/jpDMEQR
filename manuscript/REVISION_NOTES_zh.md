# 稿件说明

文章题目为 Joint Profile Inference for High Dimensional Mixed Effects Quantile Regression with Many Small Clusters。正文和补充材料分别由 `sjs_profile_quantile_main.tex` 与 `sjs_profile_quantile_supplement.tex` 生成，作者信息保存在这两个文件中。

估计参数由 unsmoothed complete ridge profile 定义。变量选择使用独立受试者，包含结局搜索与条件协变量调整。选定设计随后在两个 bandwidth 下分别完成凸优化，最终系数采用两倍半 bandwidth 估计减去完整 bandwidth 估计。方差计算保留两个估计在同一 cluster 内的协方差。

理论章节给出选择后参数与完整设计参数之间的误差界，并推导所用搜索及验证规则的统计收敛率。补充材料包含 density expansion、经验曲率界和 cluster influence representation 的证明，随机截距与时间斜率均有对应论证。

已完成的独立模拟包括24个设置，每个设置60个数据集；配对消融使用相同的1440个数据集。GSE65391计算包含52组模块工作流程，每组均完成129次受试者删除。360组独立基因分析提供1080个区间，预测部分完成900组拟合。数值表格由保存结果生成，原始结果路径见复现说明。

78个设置的完整方法比较仍待计算。既有预测比较方法目前各有一个验证折的返回记录，补充材料列出相应受试者数量和收敛状态。完整预测比较及大规模模拟的计算状态保留在复现说明中。

PDF检查读取文本与文件对象，核查引用、字体嵌入和文字页面范围。图形检查使用生成数据及程序检查，未使用视觉功能。
