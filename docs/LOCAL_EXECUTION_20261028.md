# 本地计算与保存

工作目录为 `/Users/houjian/Projects/jpDMEQR`。源码和手稿通过 Git 保存，远程仓库为 `Stork343/jpDMEQR`，计算分支为 `recovery/research-20261028`。

运行 R 时，将 `TMPDIR` 设置为工作目录中的 `results/scratch/runtime`。安装当前研究程序使用 `research/score_validity_20261026/setup.R`，依赖安装入口为 `scripts/00_install_dependencies.R`。正式实验的设置和重复次数由 `config/simulation_main.csv` 指定，方法映射及随机种子由 `config/joint_profile_submission_20261027.yml` 指定。

比较程序使用支持完整 basis 接口的 HiGHS，最低版本为1.15.1-1。依赖安装入口将其源代码版本安装到 `results/scratch/dependencies/library`，比较程序从该目录加载。独立的比较检查程序可以通过 `R_LIBS_USER` 指向同一目录。

每次重复完成后保存独立结果文件。恢复计算时，根据已经完成的重复编号继续运行，并保留所有数值失败记录。实验输出保存在 `results/raw`，汇总表和图由这些结果生成。

完整结果通过临时文件写入并改名后提供给后续程序。每个计算进程完成100次重复时，以及完成整个批次时，将完整输出复制到 `~/jpDMEQR-backups`。环境变量 `JPDMEQR_BACKUP_DIR` 可以指定其他备份位置。已完成的结果文件会保留，备份程序持续添加后续结果。

源码通过语法检查及相应数值检查后提交并推送。真实研究数据保存在本机的 `data` 目录，Git 保存程序、公开来源说明和汇总材料。

2026年10月28日的本地验证已通过66项联合推断检查和480项加权计算检查。检查结果位于 `results/preflight/recovery_20261028`。
