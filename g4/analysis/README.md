# 分析脚本与已有结果

当前目录汇集三个旧版本中的分析工具；迁移未重新生成已有 PNG，也没有把历史图作为新仿真结果。

| 文件 | 来源及用途 |
| --- | --- |
| `plot_pileup_waveform.py` | 原 g4_260623，多源总波形及已知 primary 真值分量。 |
| `compare_alpha_gamma.py`、`plot_waveforms.py` | 原 g4_260622，单粒子平均模板与尾部/总积分 PSD。 |
| `alpha_*.png`、`gamma_*.png` | 原 g4_260622 的单粒子分析图；不是多源堆积的盲重建结果。 |
| `wavelength/`、`gagg_emission_spectrum.png` | 原 g4_260621 的波长谱检查。 |

单粒子比较需要对应的数据文件；ROOT 原始数据单独保存，不在 Git 中。以下命令从仓库根目录开始，建议在单独结果目录运行，避免覆盖已有图：

```bash
mkdir -p g4/data/single-particle-check
cd g4/data/single-particle-check
python3 ../../analysis/compare_alpha_gamma.py ../gagg_waveform_20260623_15h14m49s.root ../gagg_waveform_20260623_15h14m57s.root
```

上述 γ/α 源能量分别为 661.657 keV 和 5.486 MeV。旧 ROOT 记录窗为 1 μs；平均模板绘图到 1200 ns 不代表补回了窗外尾部。两份脚本分别采用首光子对齐和不作逐事件对齐的约定，调用前应核对。单一 α 源能量也不能替代完整的 Qα(E) 标定。

本机 uproot 默认读取方式曾出现等待；需要时可显式打开 `uproot.open(path, handler=uproot.source.file.MemmapSource)` 后，将 PrimaryPhoton tree 交给模板计算函数。迁移仅修改文档中的旧目录路径，没有改变计算逻辑。
