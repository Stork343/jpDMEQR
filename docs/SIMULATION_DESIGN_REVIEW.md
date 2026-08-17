# 模拟实验设定审查报告（大规模实验前的预先准备）

审查日期：2026-08-18
审查对象：`docs/SIMULATION_PROTOCOL.md`、`config/simulation_main.csv`、`config/simulation_pilot.csv`、`manuscript/sections/main_05_simulation.tex`、`manuscript/sections/supp_06_simulation.tex`、`R/v2_simulation_v2_part*.R`、`R/v2_benchmark_adapters_v2_part*.R`、`docs/BENCHMARK_MATRIX.md`、`docs/METHOD_SPECIFICATION.md`

审查目的：在大规模实验开始前确认模拟设定充分、合理、完整、一致，避免中途返工。

---

## 〇、审阅请求（致理论部分负责人 GPT-5.6 Pro）

本报告由编码/计算团队提交，作为**模拟实验方案与理论口径的确认报告**。在冻结实验注册表、开始大规模模拟之前，请理论部分负责人逐条审阅并给出决策。每条请回答「确认 / 修改（给出修改口径）/ 裁剪（给出理由）」。

| 编号 | 待确认事项 | 需要理论方给出的口径 |
|---|---|---|
| R1（对应 P1） | 实验 E04（`lambda_gamma=4`，强 nuisance ridge）的推断目标：按 `METHOD_SPECIFICATION.md` §2，profile 目标 β\*(Λ) 依赖工作惩罚 Λ。Λ=4 时 β\* 与结构系数 β⁰ 是否显著分离？若分离，确认 E04 的 `target_mode` 由 `structural` 改为 `profile_mc`（人口逼近 β\*）；若认为仍可用 β⁰，请给出理论依据。 | β\*(Λ=4) 与 β⁰ 的分离程度判断；E04 目标口径 |
| R2（对应 P3） | 弱/强信号（a_s = 0.40 / 1.10）场景：协议 §3.2 承诺但注册表未实现。是否确认补充 2–4 个 signal 敏感性行（建议 n=200, p=1000, s=8, τ∈{0.5, 0.25}）？ | 确认补充或裁剪 |
| R3（对应 P4/P5） | copula 簇内依赖敏感性（协议 §3.4-9）与簇大小不平衡（§3.1）两个敏感性实验：是否需要做？若做，请给出 DGP 的精确设定（copula 结构/参数、不平衡分布的生成规则）。 | 做/不做 + 精确 DGP 设定 |
| R4（对应 P6） | 协议 §4.2 承诺 6 种 misspecified 场景，当前仅实现 2 种（omitted slope C06、random scale C07）。**X–随机效应相关** 与 **informative cluster size** 两个场景是否确认补充？若补充，需要理论方给出：两场景下 β\* 的明确定义（profile 目标在 misspecification 下的可识别性）、人口逼近的可行性、以及"结构系数 vs profile 目标"的评估口径。非线性随机效应贡献是否裁剪？ | 每个场景做/裁剪 + β\* 定义 + 评估口径 |
| R5（对应 S1） | Benchmark 对比矩阵：8 个适配器未实现（BIAS-ADJ-LQMM、DOUBLE-PEN-QLMM、QGEE-SCAD、QIF-SEE、TRUE-NUISANCE、POP-H、SPLIT、BAYES-MIXED-LASSO）。请确认：(a) oracle 诊断类（POP-H、TRUE-NUISANCE）必须实现；(b) 重实现型（QGEE-SCAD、DOUBLE-PEN-QLMM 等）是坚持实现，还是裁剪对比承诺并修订协议/手稿？ | 每类方法的取舍 |
| R6（对应 S3） | 人口目标逼近重复次数：协议 §7 写 2 次、手稿 main_05 写 4 次。统一为 4 次是否确认？（涉及 misspecified 场景 coverage 审计的数值误差传播） | 确认 4 次或其他 |
| R7（对应 S4） | Module A 的 p 取值：协议 {500,2000}、手稿 {500,1000,2000}、config 实际含 p=1000。确认以 p∈{500,1000,2000} 对齐三方？ | 确认对齐口径 |
| R8（对应 S2） | Module F 筛选方法：协议承诺 outcome-blind variance filter 与 no-screening 对照，注册表只有 split_quantile_score 与 CA-IQR-SIS。确认补充 variance filter 行与 no-screening 对照，还是裁剪承诺？ | 补充/裁剪 |
| R9（对应 S6/S7/S8） | 次要补强项：几何验证补 h=1.0 并强制安装 CVXR（S6）；Module G 补 target-coordinate/worker 维度（S7）；Module E 的 Λ 乘子补 0.5/2（S8）。是否全部确认？ | 逐项确认 |
| R10 | 其他：理论方对当前模拟协议（模块划分、复制数 B=1000/500、试点门槛、公平性规则、结果契约）有无需要调整之处？特别是理论章节中承诺的"有限样本理论审计诊断"（Hessian 误差、Dantzig 残差、Bahadur 余项、log-log 斜率）清单是否完整？ | 补充意见 |

> 决策原则（项目组总基准）：一切以论文的正确性、合理性、逻辑性、充分性为准，以顺利发表为目标。协议承诺过的实验场景要么兑现、要么正式裁剪并同步修订协议与手稿，不允许"协议写了但没做"。

---

## 一、总体结论

模拟实验的**骨架是完整且高质量的**：7 个模块（A–G）+ 4 个试点配置，覆盖理论尺度、误差鲁棒性、nuisance 结构、设计/精度压力、调参敏感性、超高维筛选、计算成本，主注册表 48 个实验行，数据生成机制（DGP）与协议 §3 高度吻合，调参公式与 `METHOD_SPECIFICATION.md` §9 完全一致，随机种子、失败分类、结果契约（`RESULTS_CONTRACT.md`）等防返工机制设计到位。

但是，在「**协议承诺** ↔ **注册表实际配置** ↔ **手稿描述** ↔ **代码实现**」四个层面之间发现了 **1 个会导致结论错误的实质性不一致、约 8 个协议承诺未兑现的覆盖缺口、若干文档不一致**。这些问题若在最终运行后发现，会造成大规模返工，必须在冻结注册表之前解决。

---

## 二、必须修正的问题（会导致错误结论或审稿硬伤）

### P1. E04 的 target_mode 标注错误（实质性，最优先）

- **位置**：`config/simulation_main.csv` 第 46 行（E04）。
- **问题**：E04 的 `lambda_gamma = 4`（强 nuisance ridge），notes 写明 *"Strong nuisance ridge; target changes with Lambda"*，但 `target_mode` 仍是 `structural`。
- **后果**：`METHOD_SPECIFICATION.md` §2 明确说明推断目标是**正则化 profile 参数 β\***，它**依赖工作惩罚 Λ**；Λ=4 时 β\* 与结构系数 β⁰ 明显分离。若按 `structural` 对 β⁰ 评估误差/覆盖率，会系统性出现"假偏差/假欠覆盖"，且无法与 E03（Λ=0.25）形成同目标比较。
- **修正**：E04 改为 `target_mode = profile_mc`（与 C06/C07 一致），并列入"需人口逼近"的场景清单。E03 弱 ridge 时 β\*≈β⁰ 可以保留 `structural`，但建议在手稿中注明理由。

### P2. 协议承诺的误差族/分位数覆盖未兑现（B 模块）

- **位置**：`SIMULATION_PROTOCOL.md` §6 Module B："every error family and every target quantile must appear at least twice"；"Gaussian/t3/contaminated settings at **all three quantiles** use B=1000"。
- **现状核对**（B01–B15）：
  - gaussian：0.25 / 0.75（+A 模块 0.5）✓
  - t3：0.25 / 0.5 / 0.75 ✓
  - laplace：0.25 / 0.75 ✓
  - skew_chisq3：0.25 / 0.75 ✓
  - contaminated_normal：**0.5 / 0.75，缺 0.25**（协议要求三个分位数均 B=1000）✗
  - asymmetric_laplace：0.25 / 0.75 ✓
  - **heteroskedastic gaussian：仅 0.25 一行**（不足 2 次）✗
  - **heteroskedastic t3：仅 0.75 一行**（不足 2 次）✗
- **修正**：补 3 行——contaminated_normal @ τ=0.25（B=1000）、heteroskedastic gaussian @ τ=0.5 或 0.75、heteroskedastic t3 @ τ=0.25 或 0.5。注意 B14/B15 中 heteroskedastic 的 τ 取值应覆盖 0.5 以保持分位数平衡。

### P3. 信号强度设定缺失

- **位置**：`SIMULATION_PROTOCOL.md` §3.2："weak and strong settings use `a_s = 0.40` and `1.10`"。
- **现状**：全部 48 行 + 4 行试点 `signal = 0.75`，**没有任何 0.40 / 1.10 的行**。弱信号对 TPR/FDP/功率的影响是审稿人必问的问题；强信号验证方法在乐观条件下的表现。
- **修正**：至少补 2–4 行，例如 n=200, p=1000, s=8, τ∈{0.5, 0.25}, signal∈{0.4, 1.1}（可放入 Module A 或新建一个 signal 模块行）。

### P4. copula 依赖敏感性实验缺失

- **位置**：`SIMULATION_PROTOCOL.md` §3.4 第 9 条明确要求 *"Within-cluster copula dependence: Gaussian copula with AR(1) parameter 0.4 ... used only in a dependence-misspecification sensitivity experiment"*。
- **现状**：config 无对应列/行；`r_quantile_centered_error_v2` 无 copula 实现；`generate_profile_qr_data_v2` 无 dependence 参数。
- **修正**：在 config 增加 `dependence`（或 `copula_rho`）列，DGP 增加 copula 误差生成（如 Gauss copula + AR(1) 相关结构 + 选定边际），加 1–2 个敏感性行。若决定放弃该实验，必须在协议中**删除**该承诺并说明理由——不能"协议写了但没做"。

### P5. 簇大小不平衡（imbalance）设置缺失

- **位置**：`SIMULATION_PROTOCOL.md` §3.1："An additional imbalance setting draws from a truncated shifted geometric distribution and records the realized distribution"；`supp_06_simulation.tex` 也承诺 "cluster-size imbalance" 敏感性。
- **现状**：config 无 imbalance 行；DGP 仅 `sample(m_values, replace=TRUE)` 均匀抽样。
- **修正**：DGP 增加 `m_rule`（如 `uniform` / `geometric_imbalance`）参数，config 加 1 行（如 n=200, p=1000, s=8, τ=0.5, m_rule=geometric）。

### P6. 协议 §4.2 承诺的 misspecified 场景只兑现了 2/6

- **位置**：`SIMULATION_PROTOCOL.md` §4.2。
- **现状**：协议列了 6 种 misspecified 场景：
  1. omitted random slope → C06 ✓
  2. random scale driven by random intercept → C07 ✓（`error_dist = random_scale`，代码已实现）
  3. **correlation between X and random effects → 缺失** ✗（DGP 中 X 与 b 独立生成，无相关机制）
  4. **informative cluster size → 缺失** ✗（m_i 与 X/b 无关）
  5. deliberately misspecified nuisance penalty/direction → E03/E04 部分覆盖（仅乘子，无方向错误）△
  6. **nonlinear random-effect contribution → 缺失** ✗
- **修正建议**：X–b 相关与 informative cluster size 是审稿人常规追问点，**强烈建议补上**（DGP 增加 `x_b_corr`、`informative_size` 参数 + 2 行配置，均 `target_mode = profile_mc`）。非线性随机效应贡献可裁剪并在协议注明。nuisance "方向"错误（如用错误 Λ 结构）若 E03/E04 已覆盖乘子效应，可在手稿中说明乘子即方向错误的特例。

---

## 三、建议补充/对齐的问题（充分性提升，工作量小）

### S1. Module C 比较器承诺 vs 实现状态冲突

- `SIMULATION_PROTOCOL.md` §6 Module C："Primary comparison includes LQMM, bias-adjusted LQMM, double-penalised QLMM, QGEE-SCAD, proposed and oracle proposed"。
- 但 config C01–C03 才含 LQMM/BIAS-ADJ-LQMM/DOUBLE-PEN-QLMM/QGEE-SCAD；**C04–C10 只有 PROFILE-DQR | POOLED-QR-LASSO | PROFILE-DQR-TRUE-SUPPORT**。而 `BENCHMARK_MATRIX.md` 中 BIAS-ADJ-LQMM、DOUBLE-PEN-QLMM、QGEE-SCAD、QIF-SEE、PROFILE-DQR-TRUE-NUISANCE、PROFILE-DQR-POP-H、PROFILE-DQR-SPLIT、BAYES-MIXED-LASSO 全部是 `implementation_required`（代码中返回 `benchmark_not_implemented_v2`）。
- **后果**：若这些适配器最终未实现，C01–C03 与 F01/F02 的方法列将大面积 `not_implemented`，与协议承诺的对比矩阵不符，审稿人必然质疑。
- **决策点（需项目组拍板）**：要么（a）按优先级实现适配器（推荐顺序：PROFILE-DQR-POP-H / TRUE-NUISANCE（oracle 诊断，工作量小且是理论审计关键）→ QGEE-SCAD → DOUBLE-PEN-QLMM → BIAS-ADJ-LQMM），要么（b）正式裁剪协议与手稿的对比承诺。**BENCHMARK_MATRIX.md 已规定"final run 前必须实现或显式排除"，当前配置按此规则运行会直接报错**——这正是协议想要的防返工机制，但需要现在就定实现计划。

### S2. Module F 的筛选方法覆盖不全

- 协议 §6 Module F 要求比较：无筛选（where feasible）、outcome-blind variance filter、CA-IQR-SIS、subject-split 筛选、oracle support。
- config F01–F04 只有 `split_quantile_score` 和 `ca_iqr_sis` 两列，**无 outcome-blind variance filter 行、无 no-screening 对照行、F04 也未含 PROFILE-DQR-SPLIT**。
- 建议：F04 增加一行 no-screening（p=5000 时 proposed 直接跑可能可行）；variance filter 若不做，在协议中注明"被 CA-IQR-SIS 替代并给出理由"。

### S3. 人口目标逼近的重复次数不一致（2 次 vs 4 次）

- `SIMULATION_PROTOCOL.md` §4.2/§7 说 "repeated with independent seeds"（§7 写 "independently repeated twice"）；手稿 `main_05_simulation.tex` 说 "repeated with at least four independent seeds"。
- **统一为 4 次**（更保守，符合手稿）；`approximate_profile_target_v2` 的调用脚本需实现多种子重复 + 差异检查（两次重复之差须远小于 MCSE 才可放行）。

### S4. Module A 的 p 取值三处不一致

- 协议 §6 Module A：`p in {500, 2000}`；手稿 `main_05`：`p in {500, 1000, 2000}`；config 实际：A01–A06 用 500/2000，**A07/A08 用 p=1000**（s=10）。
- 修正：以 config 为准（p=1000 的 s=10 行有存在价值），把协议 §6 Module A 的 p 集合改为 {500, 1000, 2000}，与手稿和 config 三方对齐。手稿 s∈{5,8,10} 与 config 的 s 分布（A 用 5/10、B–G 用 8/10）需在协议中注明"s 在模块间取值不同"。

### S5. C06/C07 的 p=300、C09 的 n=150 例外需文档化

- C06/C07（misspecified，需人口 MC 逼近）用 p=300：与手稿 core design 的 p∈{500,...} 不一致，但低维是合理的（population 逼近成本）。需在手稿或协议注明"misspecified 场景使用低维 p=300 以便人口逼近"。
- C09 用 n=150（簇更大 {6..12}）：手稿 n∈{100,200,400}。需注明例外（大簇时总观测数 ≈ 150×9=1350，与 n=400×5.5 相当，计算量控制合理）。

### S6. 几何验证脚本覆盖微调

- `scripts/simulation/00_validate_profile_geometry.R`：h_mult 只用 {0.75, 1.25}，建议补 1.0（协议 §11 要求 h∈{0.75,1,1.25}）；Dantzig 校验依赖 CVXR，未安装会静默跳过（只有 warning），建议在 `scripts/00_install_dependencies.R` 中强制安装 CVXR，否则几何门槛形同虚设。

### S7. Module G 计算规模维度可扩展

- 协议 §6 Module G 要求 "Vary n, p, q, selected dimension, number of target coordinates and number of worker processes"；config G01–G04 只变了 n/p/q/筛选。建议至少加 1 行改变 target coordinate 数量或 worker 数量，或注明以现有行为准。

### S8. Module E 的 Λ 乘子网格可补全

- 协议 Module E：Λ multiplier ∈ {0.25, 0.5, 1, 2, 4}；config E03/E04 只有 0.25/4 两个端点。可补 0.5/2（或注明端点足以刻画单调趋势）。

---

## 四、已核对一致、无需改动的部分（供确认）

1. **DGP 与协议 §3**：Y = Xᵀβ⁰ + Zᵀb + ε；m_i ~ Uniform{3..8}；AR(1) ρ=0.5；交替符号信号；无截距；均成立。
2. **误差分布**：8 类已实现（gaussian/t3/laplace/skew_chisq3/contaminated/asymmetric_laplace/heteroskedastic gaussian/heteroskedastic t3），全部 τ-居中（`r_quantile_centered_error_v2`，含 random_scale 用于 C07），代码正确。
3. **随机效应分布**：normal/t5/mixture/skew_lognormal 已实现且方差校正正确（t5 乘 √(3/5)、mixture 重标定、lognormal 标准化）。
4. **设计类型**：ar1（含 rho=0 即独立）/block/sparse_precision/factor/dense_precision 全部实现；dense_precision 按协议标注为压力测试。
5. **调参公式**：h=c_h·n^{−1/3}、λ_β=c_λ·√(log p/n)、μ=c_μ·(√(log p/(nh))+h²) 与 `METHOD_SPECIFICATION.md` §9 完全一致（`reference_tuning_values_v2`）。
6. **试点设计**：P01–P04（B=200）覆盖 baseline/p>n、q=2+t3、低分位+偏态、高分位+污染+随机斜率，与 §9 试点门槛匹配。
7. **复制数与 MCSE**：核心 1000 / 次要 500 / G 模块 50，符合 §7；覆盖率的 binomial 区间、配对差异的 paired MCSE 均有代码实现（`metrics_v2`）。
8. **结果契约**：run 目录不可变、失败分类 13 类、coordinate 级 schema、理论诊断 schema，`RESULTS_CONTRACT.md` 与运行框架（`_run_registry_helpers*.R`）对齐。
9. **公平性规则**：共同随机数、训练簇调参、禁止按覆盖率调参、oracle 显式标注，均写入协议与代码约束。
10. **防返工机制**：`02_run_main.R` 强制要求先跑几何验证；`BENCHMARK_MATRIX.md` 要求 final run 前实现全部 adapter（缺则 stop）——机制本身有效。

---

## 五、行动清单（建议顺序）

| 优先级 | 事项 | 涉及文件 |
|---|---|---|
| P1 | E04 target_mode 改为 profile_mc | config/simulation_main.csv |
| P2 | 补 contaminated@0.25、hetero-gaussian、hetero-t3 三行 | config/simulation_main.csv |
| P3 | 补 signal∈{0.4,1.1} 行 | config/simulation_main.csv |
| P4 | copula 敏感性：DGP 支持 + config 行 | R/v2_simulation_v2_part01.R, config |
| P5 | 簇大小不平衡：DGP m_rule + config 行 | R/v2_simulation_v2_part02.R, config |
| P6 | X–b 相关、informative cluster size 场景（或正式裁剪） | DGP + config + 协议 |
| S1 | 定 benchmark 适配器实现计划（oracle 优先） | R/v2_benchmark_adapters_v2_part02.R |
| S3 | 人口逼近统一 4 种子并写差异检查 | scripts/simulation/04_build_profile_targets.R |
| S4 | 协议 Module A 的 p 改为 {500,1000,2000} | docs/SIMULATION_PROTOCOL.md |
| S5 | C06/C07 p=300、C09 n=150 例外注明 | 手稿或协议 |
| S6 | 几何脚本补 h=1.0、强制 CVXR | scripts/simulation/00_validate_profile_geometry.R, scripts/00_install_dependencies.R |
| S2/S7/S8 | F 模块筛选对照、G 模块维度、E 模块 Λ 网格（可选） | config |

---

## 六、建议的下一步

1. 按上表先处理 **P1–P3**（纯 config 修改，半小时内可完成，可立即推送）。
2. **P4–P6 需要写 DGP 代码**（copula、imbalance、X–b 相关/informative size），改动集中在 `v2_simulation_v2_part01/02.R`，约 1 天工作量，建议与理论方（GPT-5.6 Pro）确认 misspecified 场景的 β\* 定义后实施。
3. **S1 是最大时间风险**：oracle 类适配器（POP-H、TRUE-NUISANCE）工作量小且是理论审计关键，建议先实现；QGEE-SCAD/DOUBLE-PEN-QLMM 等重实现型适配器需尽早启动或正式裁剪。
4. 所有改动落实后，再冻结 config、跑几何验证 → 试点 → 主实验。

---

*本报告基于 2026-08-18 的仓库状态（HEAD e2deced）编写。*
