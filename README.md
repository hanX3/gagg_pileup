# GAGG：g4_260621 归档版本

这是 2026-09-23 迁移时原 `g4_260621/` 工作区的快照，源码统一放在 [g4/](g4/)。原目录名作为 Git 标签保留；快照包含迁移时已有的未提交源码与小型分析结果，不表示已完整恢复旧日期的运行环境。

同一仓库还包含 [11BH_alpha_spectra/](11BH_alpha_spectra/)、[bremsstrahlung/](bremsstrahlung/) 与 [stopping_power/](stopping_power/)，这些辅助目录采用迁移时的工作区状态。原有提交经路径转换及笔记过滤后作为父提交保留，原提交号见提交说明中的 Original-Commit。

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j
cd g4/build
./gagg ../macros/run_gamma.mac
```

构建需要 Geant4 和 ROOT；旧版本是否兼容当前依赖需另行验证。输出使用相对路径 `../data`，因此从 `g4/build/` 运行。历史物理参数以此版本源码为准，旧数据运行条件以对应 ROOT 的 RunInfo 为准。

ROOT 原始数据、构建产物、AGENTS.md、jupyter/ 和所有 .ipynb 不纳入 Git；本地完整原目录另行归档。使用 `git diff g4_260622 g4_260623 -- g4/` 比较版本，当前开发和完整说明见 main 分支。
