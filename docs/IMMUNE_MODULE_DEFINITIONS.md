# Immune module定义

五个module沿用G2血液转录组框架的既定功能分类。成员与计量单位已经在查看expression与C3关联以前确定。`config/module_signatures.csv`保存五个module的完整gene表示，共180行；`config/module_score_members_20261022.csv`保存本项目评分使用的完整成员，共181行。两个文件各有独立的`.freeze.csv`，评分规则位于`config/module_scoring_20261022.yml`。

## 公开来源

[Banchereau et al 2016 Cell](https://pmc.ncbi.nlm.nih.gov/articles/PMC5426482/)使用已有260个血液module分析GSE65391，并引用[Chaussabel and Baldwin 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4118927/)。后者的Dataset analyses明确说明，可以对全部module成员的normalized expression计算平均值。Banchereau文章的Modular transcriptional analysis以probe为统计单元计算百分比活动分数，其Figure 7也使用immune signature的平均normalized expression。本项目将完整成员的等权表达均值作为quantile regression的预定covariate。

原始G2 probe成员来自[Hejblum et al 2015的S7公开附件](https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1004310.s007)。附件README将`GMTs_PLOScb.RData`中的`gmt_modulesV2`明确指向原G2作者wiki，RData保存260个module的ILMN probe列表和功能名称。该文件日期为2014年10月20日。完整gene表示来自公开维护的[tmod 0.50.13 source archive](https://cran.r-project.org/src/contrib/tmod_0.50.13.tar.gz)，其`tmod.rda`以`SourceID=DC`、`original.ID=M*`记录G2 module，并保留原作者wiki的module链接。[tmod说明文档](https://www.stats.bris.ac.uk/R/web/packages/tmod/vignettes/tmod.html)解释了DC集合的来源和gene表示。

完整文章、Cell补充PDF、S7 ZIP、tmod archive以及对应解析文本均保存在`data/raw/public_modules/banchereau2016_20261022/`。该目录的`source_manifest.csv`列出10份已完整读取或解析的来源文件、公开URL、文件大小与保存时间。模块定义使用G2既定功能名称，包含IFN的M1.2、M3.4、M5.12，Plasma Cells的M4.11及Neutrophils的M5.15。

## 官方平台注释

平台来源为[GPL10558官方annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz)及[原始submitter表](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL10558&form=text&view=full)。两份完整来源均由GEOquery解析，各有48107个probe，ID集合一致。curated annotation日期为2016年8月9日。原始submitter的Species字段用于区分human transcript与array control，Symbol、Entrez_Gene_ID及Synonyms字段保留在注释输出中。

`data/interim/GSE65391/probe_annotation_20261022.csv`覆盖expression全部43799个probe。计数为41181个unique_gene、13个multiple_genes和2605个gene_unannotated；这些观测probe均属于human transcript。gene评分保留unique_gene；M5.15直接使用既定的24个ILMN probe测量。既有expression、metadata与空注释文件保留原状。

三个源LOC名称可以通过NCBI Gene ID获得唯一的正式名称，且GPL10558存在相同Gene ID及正式symbol的观测probe。

| 源成员 | NCBI Gene ID | 正式symbol | 观测probe示例 |
|---|---:|---|---|
| LOC129607 | [129607](https://www.ncbi.nlm.nih.gov/gene/129607) | CMPK2 | ILMN_1783621 |
| LOC26010 | [26010](https://www.ncbi.nlm.nih.gov/gene/26010) | SPATS2L | ILMN_1683678 |
| LOC728358 | [728358](https://www.ncbi.nlm.nih.gov/gene/728358) | DEFA1B | ILMN_1679357 |

官方记录的完整JSON保存在`ncbi_gene_summary.json`，包含human taxonomy、current official nomenclature及Gene ID。本地`public_member_mapping.csv`逐项保留source_gene_symbol、gene_symbol、source_gene_id、source_mapping_url与mapping_status。180个源gene成员全部具有唯一观测映射，源成员数量保持完整。

## 成员覆盖与固定表示

下表的probe计数直接使用`expression_processed_probe.rds`的rownames核对。gene与probe分别计量。

| module_id | 功能名称 | 完整gene表示 | 原始probe覆盖 | 评分表示 | 固定权重 |
|---|---|---:|---:|---|---:|
| M1.2 | IFN | 27/27 | 30/36 | 27 genes | 1/27 |
| M3.4 | IFN | 53/53 | 54/63 | 53 genes | 1/53 |
| M4.11 | Plasma Cells | 18/18 | 17/20 | 18 genes | 1/18 |
| M5.12 | IFN | 59/59 | 50/64 | 59 genes | 1/59 |
| M5.15 | Neutrophils | 23/23 | 24/24 | 24 ILMN probes | 1/24 |

每个module使用全部固定成员，`minimum_fraction=1.0`。源gene与当前symbol的唯一名称转换具有记录，原始probe成员按原ID保存。M5.15的score保留原始24个probe及相同权重。该module的23个gene表示也保存在来源表中，评分文件明确标记`member_type=ILMN_probe`。

## 评分计算

gene表示使用processed Series Matrix中已经normalized的表达值。多个probe映射至同一gene时，按training samples的MAD选取一个probe，既定次序下的首个probe处理MAD相同的情况。每个module计算完整成员的等权均值，随后使用training rows计算该module均值的center与sample SD，转换全部samples的分数。M5.15直接计算24个既定probe的等权均值，然后进行相同的training转换。

固定成员缺少测量、表达出现非有限值、权重改变或training SD为零时，程序立即终止。输出保留每个module的计量单位、成员数量、完整权重、gene到probe的选择，以及training center和SD。评分过程仅读取expression及预定training记录。

`R/v2_application_v2_part06.R`实现来源冻结检查、官方注释合并和固定表示评分；`prepare_gse65391_lowdim_design_v2`接收独立的`module_score_members`。五个score与clinical covariates组成的设计在129名患者、867次就诊中具有完整rank，Model A为16维，Model B为17维，既定上限为26维。

## 全流程sample standardisation

`R/v2_application_v2_part08.R`将preprocessing纳入低维paired估计器。对给定subject样本，重新计算clinical时间center、缺失值处理和类别编码，使用该样本选择各gene的最高MAD probe，再为全部五个module计算固定成员均值及其sample center和SD。module coefficient使用该样本一个module SD的单位。固定源成员、module计量单位与成员权重保持既定定义。

每次leave-one删除完整subject，并在余下样本重新调用相同的prepare函数，然后计算paired one-step。leave-one coefficient使用该次样本的standardisation，coordinate名称与完整样本一致。jackknife将这个sample-standardised估计器作为整体计算，全部拟合使用完整样本的固定h及h/2。每次保存clinical处理、gene到probe的选择、全部module center和SD、module audit及完整fit对象。

输出字段`jackknife_scope=complete_sample_standardised_workflow`标记这个预定敏感性区间，`sandwich_scope=profile_score_given_full_sample_preprocessing`标记完整样本给定design下的paired sandwich。其target保留`unsmoothed_regularised_profile`及sample-standardised parameterisation。数据自适应probe选择和standardisation的研究条件由application-matched模拟与相应理论检查处理。

## 可重复检查

`06_build_immune_modules.R audit-expression-probes`生成逐probe观测表与五模块覆盖计数。`check`从S7、tmod和完整官方GPL来源重新建立gene定义，核对当前冻结文件。`freeze-scores`保存固定评分成员及其冻结记录。运行`07_check_immune_modules.R`可检查完整表达均值、五模块覆盖以及clinical/module设计rank。

`test-application-module-sources.R`包含20个通过的assertions。测试使用真实source文件和expression，逐项核对完整成员、LOC恢复与直接表达平均值，并以相同training记录、不同validation记录检验training转换保持一致。来源解析使用GEOquery、R的RData读取与tar/ZIP工具；HTML和xlsx检查使用BeautifulSoup与openpyxl，Python依赖位于`scripts/application/requirements-public-modules.txt`。

`08_check_lowdim_workflow_GSE65391.R`声明两个固定坐标用于研究计算检查，按两个source batch各选择六名真实subject，逐次重新计算全部五个module的preprocessing。`test-application-profile-workflow.R`的63个assertions通过，涵盖完整subject删除、相同bandwidth与coordinate名称、全部module audit、每次独立training记录，以及直接重新准备设计和拟合。完整应用入口使用全部既定坐标。
