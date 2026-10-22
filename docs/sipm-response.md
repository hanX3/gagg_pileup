# S13360-6050PE 光子到电信号响应

当前统一采用 **Hamamatsu S13360-6050PE，25°C、3 V 过压**，配置为 [sipm_hamamatsu_s13360_6050pe_v2.json](../g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)。原生宏、批量运行器、离线验证、器件分析和当前笔记均使用这一器件。参数状态 `datasheet_based` 表示厂家典型值及注明的模型近似；实物标定尚未完成。

文献关系见 [L06–L07](literature-use.md#sipm-use-20261010)，器件来源与采用边界见 [L13](literature-use.md#hamamatsu-s13360-20261013)，完整公式、Markdown参数表与波形见 [SiPM笔记](../jupyter/codex_0001_SiPM_Digitization.ipynb)。

<a id="native-sipm"></a>
## Geant4 原生数字化

通过 `G4VDigitizerModule` / `G4DigiManager`，在完整窗的光子按到达时间排序后显式数字化，保存 `SiPMEvent` 到同一 ROOT；保留原有光子记录。输入为 t/x/y/光子能量，α/γ标签不参与器件响应计算。

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j 2
ctest --test-dir g4/build --output-on-failure
cd g4/build
./gagg ../macros/run_single_alpha_sipm.mac
./gagg ../macros/run_single_gamma_sipm.mac
```

构建优先读取 nlohmann_json CMake 包；若仅安装了带 `nlohmann/json.hpp` 的ROOT，可直接使用该头文件并核查版本≥3.2。单粒子宏均启用S13360-6050PE v2，原生模式要求 `/run/beamOn 1`。20 μs窗、5 μs发射、中央[5,15) μs、中央正入射；α源能5.486 MeV，γ源能661.657 keV。

### 器件与读出设置

| 项目 | 当前值与依据 |
| --- | --- |
| 有效面 / pixel | 中央6 × 6 mm；120 × 120＝14400个50 μm pixel。 |
| 温度 / 过压 | 25°C / 3 V；厂家典型工作条件。 |
| 满恢复增益 / 电荷 | 1.7 × 10⁶ electrons；q₀=0.272370 pC。 |
| PDE | 450 nm约40%，530 nm约35%；近似读取PE曲线，已含填充因子。 |
| 暗计数 / 串扰 | 2 MHz / 3%；分别作初级Poisson率和简化单邻居递归概率。 |
| 恢复 | 29 ns，同系列CS器件代理，尚非6050PE实测值。 |
| AP / SPTR | 暂不附加，零参数表示未建模。 |
| 电流核 / 跨阻 / 低通 | 2/20 ns双指数差、50 Ω、2 ns一阶低通，待实际读出标定。 |
| 采样 / 孔径 | 2 ns / 0.5 ns；内部步长0.25 ns。 |
| 基线 / 噪声 | 10 mV / 0.1 mV RMS。 |
| ADC | 12 bit，0–1000 mV，LSB=0.244140625 mV。 |

来源为已保存的[滨松KAPD1052E07数据表](https://www.hamamatsu.com.cn/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360_series_kapd1052e.pdf)，PDF2–5页。PDE限定400–650 nm、0–3 V：3 V近似读图，2 V按450 nm处0.315/0.400缩放，更低过压线性连零。波长/过压可分离与低过压恢复属于项目近似，说明性读图精度约±1个百分点；不外推、不重复乘74%填充因子。

29 ns来自Rosado2018 PDF3表1的S13360-1350CS/3050CS，不能写成6050PE厂家指标。真实封装、光学脂、AP、SPTR和SPE/前端仍待标定。10×10 mm虚拟层保存完整到面光子，数字化只接受6×6 mm有效面；该记录层并不代表器件尺寸。

### 数据与随机数

`SiPMDigitizer`生成G4Event所有的digi集合；`SiPMResponse`是项目数值核心。原生首版只支持`double_exponential`；下文`measured_voltage`为Python接口。输入光子保持稳定时间顺序，整窗共享pixel恢复历史。

雪崩使用事件种子派生的5条独立`std::mt19937_64`流，电子学噪声另用读出种子，不消耗CLHEP输运随机数。跨编译器/NumPy不要求逐次随机输出相同；对同一雪崩表的无噪读出可逐点比较。

| 树 / 字段组 | 内容 |
| --- | --- |
| RunInfo | 完整`sipm_config_json`、模块版本、基础种子和输运/窗口条件。 |
| WaveformEvent / PrimaryPhoton | 完整到面光子和粒子真值，保持原始记录。 |
| SiPMEvent：`avalanche_*` | 时间ns、微胞、类型、父编号、原光子索引、恢复系数、电荷C。 |
| `current_time_ns/current_A` | 解析格内平均电流，尚未过读出低通。 |
| `sample_time_ns/signal_mV/analog_mV/adc_counts/adc_clipped` | 孔径中心、无噪电压、含基线噪声电压、ADC与削顶。 |
| `charge_total_C/charge_gate_C/current_integral_C/current_tail_C` | 全窗与按触发时刻选门的雪崩电荷、窗内电流积分、窗后解析尾部。 |
| `n_*`、`avalanche_seed/readout_seed` | 机制计数、削顶/窗外候选计数、实际事件派生种子。 |

### 当前器件的已保存结果

500 mV诊断6窗、1000 mV复跑5窗保留在`jupyter/run_records/sipm_native_20261013_095254/`和`sipm_native_20261013_095833/`。α源能5.486 MeV、沉积4.939493 MeV，38453个记录面光子中31308个进入有效面；条件γ源能661.657 keV、沉积0.202202 MeV，7370个光子中5966个进入有效面。无沉积γ也完整保留。

α/γ模拟峰值508.018/180.496 mV；500 mV时α有4点削顶，1000 mV配置无削顶，同种子量程对照的光子、雪崩及模拟波形一致。α有效面满恢复加权PDE为34.943%，光子触发10367次。恢复20/29/40 ns单例敏感性为+2.47%/0/−2.76%，不表示实物误差范围。

当前软件检查包括14项C++基准、31项Python检查、逐ROOT审计、随机数隔离、C++/Python读出与内部步长收敛，以及14400 pixel短闪光的解析占用参考。它们验证模型实现，不能代替器件标定或给出能量分辨率/堆积效率。

固定样本与三组随机数对照从项目根目录运行：

```bash
.venv/bin/python -m analysis.single_particle_sipm run --alpha 4 --gamma 32
.venv/bin/python -m analysis.single_particle_sipm analyze jupyter/run_records/<本次目录>
```

运行器使用预定种子、保留全部无沉积γ；单γ宏使用已知有信号种子261010518，属于条件展示。源码、实际配置/宏、随机种子、环境、Git改动与输入输出哈希均保存。每ROOT一窗，串行启动间隔至少1秒，ROOT写入`g4/data/`，记录与图分别写入`jupyter/run_records/`和`jupyter/artifacts/`。

只分析指定器件的PDE、波形、占用和恢复敏感性：

```bash
.venv/bin/python -m analysis.hamamatsu_s13360 \
  jupyter/run_records/sipm_native_20261013_095833 \
  --output jupyter/artifacts/s13360_characterization_new
```

输出目录必须是新目录；模块读取该批实际配置和ROOT，不依赖其他器件的历史基准。

### α / γ 归一化波形比较

[sipm_waveform_comparison.py](../analysis/sipm_waveform_comparison.py) 比较同一器件配置下的两个单粒子 ROOT，分别叠加有效面光子计数、SiPM 电流、analog 和 ADC 电压。电压先扣除 RunInfo 记录的基线，各阶段独立除以 0–2000 ns 内的峰值；各事件以首个有效面光子到达为零，同一事件四阶段共用零点。原始光子用 2 ns bin，保留电信号中的暗计数与噪声。

取得对应 ROOT 后，在项目根目录的 Python 环境中运行：

```python
from analysis.single_particle_sipm import load_root
from analysis.sipm_waveform_comparison import plot_normalized_comparison

examples = {
    "alpha": load_root("g4/data/20261013_09h58m53s.root"),
    "gamma": load_root("g4/data/20261013_09h59m09s.root"),
}
summary = plot_normalized_comparison(examples, "jupyter/artifacts/sipm_normalized_example")
```

输出四阶段 PNG 与数值 JSON；重复调用更新指定目录下的这两份派生文件。当前这两个样例中，α 较长的尾部在 SiPM 前后均保留，analog 与 ADC 几乎重合。摘要同时给出 100–2000 ns 相对 0–2000 ns 的积分占比；该时间门仅用于描述形状，未经 PSD 优化。两例的沉积能量、光子统计和微胞占用不同，峰值归一化不会消除这些差异；首光子对齐使用仿真输入，不能据此给出总体粒子甄别性能。图表与完整说明保存在本地 SiPM 笔记第 5 节。

## 离线运行

工作目录为项目根目录，使用已有 `.venv`；不增加第三方依赖。

```bash
.venv/bin/python -m unittest pileup.analysis.test_sipm_response -v

.venv/bin/python -m pileup.analysis.sipm_response \
  g4/data/20261010_09h47m20s.root \
  --config g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json \
  --output jupyter/artifacts/sipm_example_001 \
  --seed 2610101011 --readout-seed 2610101012
```

ROOT 是本地字段回归样本，不随 Git 分发。替换成自己的新 ROOT 即可；指定输出目录必须尚不存在。所有参数和路径由命令明确传入，不批量改写已有光子数据或冻结实验。

完整的合成检查、三种响应设置比较和读出时间步收敛：

```bash
.venv/bin/python -m pileup.analysis.validate_sipm_response \
  g4/data/20261010_09h47m20s.root \
  --config g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json \
  --output jupyter/artifacts/sipm_validation_001
```

验证运行器读取 S13360-6050PE 配置，生成受控双闪光恢复、瞬时 Poisson 照明占用基准，以及同一 ROOT 的 `linear`、`recovery`、`full` 三组输出。`full` 保持器件配置；前两组只关闭注明的效应。解析极限检查分别临时采用3000个胞/PDE=1与256个胞/PDE=0.35，关闭暗计数、相关噪声和渡越延迟；这些条件写入结果说明，不当成真实器件预测。还从真值树选择一个距窗末足够远的 α，将其光子历史首光子对齐，构造单 α 及间隔0/30/100/300/1000 ns的双拷贝输入；保存所选真值、输入光子、实际配置及结果。双拷贝共享同一光学历史，不是两次独立 Geant4 事件；这一已知构造只检验响应，不是盲分离或器件标定。

## 输入、概率与时间顺序

输入要求一个 ROOT 含一个 `WaveformEvent` 和一个 `RunInfo`。读取已排序且索引对应的 `all_photon_arrival_t_ps/x_mm/y_mm` 与 `photon_energy_eV`，转换到 ns/mm/eV。时间必须在 `[0,T)`，T 与中央区从当前文件的 RunInfo 读取。

- 支持一个读出通道、矩形规则微胞布局。x/y 映射到完整 pitch 单元，器件 PDE 已含的微胞填充因子不重复扣除。器件外的光子计入 `outside_footprint`。外轮廓最大坐标归入最后一列/行，内部边界归入右侧/上侧单元。
- `pde_mode=table` 使用二维 `pde_probabilities[过压][波长]`，做双线性插值；表需覆盖零过压到工作过压以及有效入射光谱，零过压行为零，禁止外推。波长由 `1239.8419843320026 / photon_energy_eV` 得到 nm。
- `pde_mode=constant` 的 `constant_pde` 是满恢复有效 PDE；恢复开启时再乘配置的相对触发曲线。旧 ROOT 缺少能量时还需明确传 `--allow-missing-energy`，并在清单标记该近似。不能对缺少能量的 ROOT 使用二维 PDE 表。
- 初始全部微胞充分恢复，无窗前后脉冲。成功雪崩后按 `r = 1-exp(-Δt/recovery_ns)` 恢复，电荷为 `q0*r`，`q0 = gain_electrons*e`。失败候选不重置状态。`recovery_enabled=false` 是去除占用/恢复影响的线性对照。
- 对光子，二维 PDE 表直接提供当前波长/过压下的触发概率，只抽样一次。光子能量不与单次雪崩电荷成正比。
- 常数 PDE 模式和非光子候选使用相对曲线 `f(r) = [1-exp(-max(r*U-U0,0)/Uch)] / [1-exp(-(U-U0)/Uch)]`。`U=overvoltage_V`，`U0=trigger_threshold_V`，`Uch=trigger_scale_V`。它是明确的经验子模型；噪声载流子的触发曲线不假定由光谱 PDE 表自动推导。

光子以有序数组读取，暗事件与二次事件使用最小时间堆，共用 `last_fire_ns`。同一时刻，已入队的暗/二次候选先于尚未处理的光子；队列按入队顺序打破并列，光子按输入稳定顺序打破并列。该约定在1 ps输入精度下确定性重现，不代表已经模拟亚 ps 内的微观竞争。

## 暗计数、串扰、后脉冲及 SPTR

所有机制可分别关闭。每种雪崩都更新微胞，并可产生下一代二次候选。

| 参数 | 本实现的精确定义 |
| --- | --- |
| `dark_rate_hz` | 整颗器件在充分恢复时的初级暗雪崩率，不含串扰/后脉冲。按 Poisson 过程生成候选并均匀分配微胞，目标接受概率为 f(r)。受照明占用时实际暗雪崩率可降低。 |
| `crosstalk_probability` | 一次满电荷雪崩产生一个串扰候选的概率。部分电荷时乘母雪崩 q/q0；从现有的上下左右相邻微胞均匀选目标，再按目标 f(r) 接受。每次母雪崩至多一个候选，但允许递归。 |
| `crosstalk_delay_ns` | 为零时 prompt；正值时为候选延迟的指数分布均值。 |
| `afterpulse_probability` | 一次满电荷雪崩产生一个同胞陷阱释放候选的概率。部分电荷时乘 q/q0，发生时按目标 f(r) 接受。不是已经包含恢复损失的实测后脉冲计数比。 |
| `afterpulse_delay_ns` | 单指数陷阱释放延迟均值。每个母雪崩至多一个候选，允许其成功后继续产生下一代。 |
| `transit_delay_mean_ns/sigma_ns` | 外来光子到触发候选的因果延迟模型。sigma>0 时用正值 Gamma 分布，shape=(mean/sigma)²、scale=sigma²/mean；sigma=0 时用固定延迟。sigma 是此假设分布的标准差，不是 FWHM。延迟在重新排序及竞争之前加入。 |

这里的串扰/AP 是便于检验的单候选分支模型，不声称涵盖真实器件的空间核、多陷阱时间常数或全部噪声机制。测得的多光电子比例、平均额外电荷 c 或含相关噪声的阈值 DCR 不能不经拟合直接填入上述概率。要求串扰与 AP 概率之和小于1，保证候选分支的上界为次临界；超出 `max_candidates` 时失败并记录原因，不返回静默截断的有效波形。

温度字段记录参数的适用条件，程序不会自动由温度推算 PDE、增益、DCR。改变工作过压时也需同步提供对应增益和噪声参数。不同机制使用从雪崩种子派生的独立随机流；读出噪声另用 `--readout-seed`。不保证不同配置具有逐事件一一相同的随机路径。

## 单雪崩响应与 ADC

`double_exponential` 使用面积归一电流核：

```text
k(t) = [exp(-t/decay_ns)-exp(-t/rise_ns)] / (decay_ns-rise_ns), t>=0
I(t) = sum(q_k * k(t-t_k))                 # 换算时间后单位 A
V(t) = I(t) * transimpedance_ohm           # 可接一阶 RC 低通
```

`rise_ns < decay_ns`，输出衰减时间和微胞恢复时间独立配置。该核的各时间格电荷积分根据原始连续雪崩时刻解析计算，不先四舍五入雪崩时间。结果 `current_A` 是格内平均电流。`lowpass_ns>0` 时，以格内常值电流近似驱动一阶 RC，记录零初始状态；须检查 `internal_dt_ns` 收敛。核依据 Acerbi 的形式，但面积归一和这里的电路为项目建模选择。

Python入口的 `measured_voltage` 读取已去基线、对应满恢复单雪崩、**包含实际读出链**的 CSV：

```csv
time_ns,voltage_mV
0,0
...
```

时间从0起严格增加。模板可包含负向或双极性部分；范围之外视为零，不补猜长尾。程序按 q/q0 缩放并叠加，忽略 `rise_ns/decay_ns/transimpedance_ohm`，要求 `lowpass_ns=0` 以避免重复同一电子学。CSV 路径相对于配置文件目录；运行时复制模板及其哈希。模板模式中雪崩时刻向后对齐到内部网格，附加延迟小于一个格宽，需做网格收敛。单雪崩模板已含时间抖动时，应单独评估是否又启用了相同 SPTR。

ADC 在每个 `sample_period_ns` 的起始处对 `aperture_ns` 范围平均，时间标签为孔径中心。周期、孔径及窗口需与内部网格整除，孔径不超过周期。然后加 `baseline_mV` 和采样后的独立高斯电压噪声 `noise_rms_mV`，按 `floor((V-Vmin)/LSB)` 编码并削顶，`LSB=(Vmax-Vmin)/2**bits`。这是明确的采样孔径和噪声近似，尚无实测有色噪声、AC耦合、放大器大信号非线性或抗混叠验证。

## 输出和复现

| 文件 | 内容 |
| --- | --- |
| `avalanches.npz` | `records` 结构化数组，以及各微胞末次成功雪崩时间 `last_fire_ns`。未触发微胞用 `-inf`。 |
| `waveform.npz` | ADC采样时刻 `time_ns`，无噪信号 `signal_mV`，含基线噪声的 `analog_mV`，`adc_counts`，削顶标记 `saturated`；解析电流模式另有 `current_time_ns/current_A`。 |
| `manifest.json` | 配置、种子、输入ROOT哈希/RunInfo、代码版本、源码哈希、执行环境、计数、窗外候选、积分及电荷尾部、耗时和完成/失败状态。 |
| `config.json`、`source/` | 实际配置和数字化/读出/测试源码副本；实测模板另复制为 `spe_template.csv`。 |
| `git_status.txt`、`source_changes.patch` | 工作区改动状态和已跟踪数字化源码的未提交差异；尚未跟踪源码的完整内容由 `source/` 保存。 |

雪崩数组字段：`time_ns, cell_id, charge_C, kind, parent_id, photon_index, recovery_fraction`。数组行号即雪崩编号；`kind` 0/1/2/3 分别表示外来光子、初级暗事件、串扰、后脉冲。`parent_id=-1` 表示无母雪崩，`photon_index=-1` 表示该事件不是直接由 ROOT 中某个外来光子触发；沿父编号可追溯相关事件。光子索引总是输入完整窗中的原索引，包括引入触发延迟重新排序之后。

雪崩次数、电荷/q0 和 ADC 积分分别保存与解释，部分恢复时它们不相等。manifest 中 gate_charge 是**触发时刻位于中央区的雪崩电荷和**，不是中央区内电压或电流积分，后者还包含跨区尾部。

从完整窗口起点推进并保留保护区，不在中央区起点重置状态。独立 ROOT 分别初始化，不串接成连续采集。窗末之后的二次候选按类别计数而不继续展开；解析电流记录未落入窗口的尾部电荷。现有5 μs前保护区是否足以消除初态影响，必须结合实际器件和光学输入另做收敛检查。

采用真实 SiPM/电子学后需重新建立响应模板及电荷—沉积能量标定。当前没有改变 G4 中的虚拟吸收面；封装反射、返回晶体的光路等仍需要对应光学模型和新的 Geant4 运行，离线 PDE 不能补出这些轨迹。
