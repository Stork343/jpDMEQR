# SLE应用计算设计

真实SLE输入保持129名subject、867次就诊，结果为观测log C3。固定module为M1.2、M3.4、M4.11、M5.12、M5.15，来源与成员见`IMMUNE_MODULE_DEFINITIONS.md`。模拟结果使用独立目录和`declared_matched_simulation`标记，原始C3及实际数据文件保持原值。

## 计算矩阵

Model A和Model B在tau=.25、.50、.75进行完整clinical/module推断，共六个规格。每个规格计算完整样本及129次全流程leave-one，合计780次paired workflow拟合。Lambda固定为1。Lambda=.5和2分别形成六个规格的敏感性计算，h multiplier=.75和1.25也按相同方式计算；random-intercept加time-slope使用六个规格。tau=.10仅使用曲率检查通过的h multiplier=1和1.25。各项计算保留全部五个module。

module报告pointwise 95%区间，以及覆盖两个model、三个既定tau和五个module的Bonferroni family30区间。clinical系数独立保存。全部jackknife区间记录为完整sample-standardised workflow的敏感性计算，各次重新准备设计后使用完整样本的固定h及h/2。

tau=.10使用独立的两个model与五个module组成family10。不同h、Lambda与random-effect设计分别保存family身份。

repeated validation使用五次重复、每次五个outer folds。每个outer training set内包含五个subject folds，全部preprocessing在相应training subjects中完成。profile预测使用两遍score loadings的初始convex L1系数，h_est multiplier的三个数值、Lambda的三个数值、lambda multiplier的四个数值与screen dimension的三个数值构成108个候选。候选通过equal-subject inner pinball loss选择，数值相同时按固定candidate顺序处理。filter sizes为2000、5000、10000；random-slope敏感性使用5000。

探索层使用十个预定50/50 subject splits并交换角色。selection侧选取三个gene，独立inference侧使用clinical调整与这三个gene。逐split保存目标和区间，汇总选择频率及系数median和IQR。

## Matched总体与数据生成

`application_matched_prospective_20261022.yml`预先给出两个model、q1/q2、normal/t3和三个tau的24个情景，每个情景B=200。真实129个subject的经验分布定义固定总体preprocessing。总体gene的probe选择、module center和SD均固定；生成结果使用该固定X及声明的系数。每次重复重新抽取129个完整subject，并重新计算全部preprocessing。

intercept取log100。time、age与SLEDAI按总体SD的系数分别为−.10、.05及−.05；五个module的系数依序为−.15、0、.10、.05及−.10。其余categorical系数为零。normal和t3误差的边际SD均为.20，t3除以sqrt3后缩放。normal cluster intercept的SD为.15，q2的slope按总体time SD缩放至.15，intercept与slope独立。这些参数定义模拟研究的生成过程。

tau=.50时，连续中心对称误差和固定ridge penalty使nuisance minimiser与profile score具有奇对称性，声明的生成系数为profile目标。tau=.25使用100000个cluster的四次独立完整profile计算，保存score、curvature、nuisance、MCSE及half-bandwidth refinement。tau=.75通过精确反射`2*beta0-beta_star(.25)`获得，保留derived身份，并在独立总体样本计算score audit。

全部模拟系数转换到固定总体单位。module斜率乘以`SD_population/SD_replication`，intercept加入module center与clinical time center的调整。probe选择变化会改变predictor测量，记录其数量和measurement scope；此时转换提供共同单位，流程差异由实际模拟表现呈现。每次leave-one也重新进行单位转换，然后计算整体jackknife变化。

## Comparator定义

POOLED-QR-LASSO调用`quantreg::rq.fit.lasso`，使用官方`LassoLambdaHat`的training plugin及既定四个倍率。clinical/module的penalty为零。QGEE-SCAD复用geeVerse0.3.1的官方`.fit_single_lambda` kernel，按原30点lambda路径和返回HBIC选择，保存全部path、converged和iterations。其输入按subject与time排序，按官方intercept定义提供penalty factor。

LQMM调用官方lqmm包，在完整clinical/module设计上使用七个quadrature points和既定likelihood tolerance。公式中的临时列名与原coordinate分别保存。三种comparator均使用自身的统计目标，profile几何诊断字段标记为不适用。真实输入的source-specific对照检查在拟合前完成。

## 构建身份与输出

subject folds、inner folds和探索角色保存于独立derived目录。完整关联与validation读取实际matching method gate，保存构建commit。目标asset分别保留自身的构建commit；承接目标时核对同一method version、完整DGP config以及固定总体X、Z、subject rows和preprocessing，当前拟合构建身份另行记录。

population construction与数值模拟使用独立核心检查。程序读取实际compiled 112项与population 319项检查结果，要求failure、error、skip及warning均为零，并核查安装DLL的两个注册入口、已提交的matched定义和完整总体设计。这项检查的scope明确为population construction及simulation，记录应用构建与核心构建commit。实际C3关联仍读取完整matching application gate。

每个matched repetition保存全流程checkpoint、总体单位下的系数和两种标准误，以及target MCSE和区间覆盖。CV保存全部候选的subject loss、选择结果和outer predictions。所有应用表格从这些机器可读结果生成。

应用source-specific检查共176个assertions通过，包含真实数据的profile计算、完整subject jackknife、source成员和评分、DGP来源保留与单位转换，以及三个comparator的官方source调用对照。
