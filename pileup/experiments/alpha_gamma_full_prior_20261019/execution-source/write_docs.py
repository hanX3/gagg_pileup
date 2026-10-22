from pathlib import Path
import json

project=Path(__file__).resolve().parents[3]
folder=project/'pileup/experiments/alpha_gamma_full_prior_20261019'
artifact=project/'jupyter/artifacts/alpha_gamma_full_prior_20261019'
dev=json.loads((artifact/'development/summary.json').read_text())
test=json.loads((artifact/'test/summary.json').read_text())
diag=json.loads((artifact/'test/transition-diagnostics.json').read_text())
labels={'baseline':'扩展标定物理拟合','unet':'U-Net','hybrid_none':'混合：无显式惩罚',
        'hybrid_both':'旧混合：Q＋类别','class_only':'新混合：仅类别'}
def lookup(summary,method,particle):return next(r for r in summary['pooled'] if r['method']==method and r['particle']==particle)
def recovery(row):return f"{row['recovered_count']}/{row['target_count']}（{100*row['full_recovery']:.1f}%）"
def table(summary):
    lines=['| 方法 | α完整恢复 | γ完整恢复 | α完整precision | γ完整precision | 额外预测 |',
           '| --- | ---: | ---: | ---: | ---: | ---: |']
    for method in summary['methods']:
        a,g=lookup(summary,method,'alpha'),lookup(summary,method,'gamma')
        lines.append(f"| {labels[method]} | {recovery(a)} | {recovery(g)} | {100*a['full_precision']:.1f}% | {100*g['full_precision']:.1f}% | {summary['diagnostics'][method]['unmatched_predictions']} |")
    return '\n'.join(lines)
failed=[r['check'] for r in test['frozen_screen']['checks'] if not r['passed']]
decision=('本批独立数据通过与开发相同的筛选。仍只报告本批有限窗口的证据，暂不替换默认拟合流程。' if not failed else
          '本批独立数据未通过与开发相同的筛选，不替换默认拟合流程。未通过项：'+', '.join(failed)+'。')
contrast='\n'.join(f"- 相对{labels[r['b']]}的{r['particle']}恢复数变化{r['delta_count']:+d}，恢复率差{100*r['delta_recovery']:+.2f}个百分点；成对整窗95%范围[{100*r['paired_window_p025_p50_p975'][0]:+.2f}, {100*r['paired_window_p025_p50_p975'][2]:+.2f}]个百分点。" for r in test['paired_contrasts'] if r['b']=='baseline')
transition=next(r for r in diag['transitions'] if r['reference']=='baseline' and r['particle']=='gamma')
weak=next(r for r in diag['strata'] if r['group']=='weak_gamma_500_1000')
weaktext='、'.join(f"{labels[m]} {weak['outcomes'][m].get('recovered',0)}/{weak['n']}" for m in ('baseline','hybrid_both','class_only'))
body=f'''# α＋γ完整混合流程：拆开显式电荷/类别惩罚（2026-10-19）

## 研究问题与唯一改动

前轮[固定候选消融](../alpha_gamma_prior_ablation_20261019/README.md)提示跨类型共用网络Q惩罚可能损害γ。本轮恢复原完整流程的时刻搜索与残差补峰，仅将显式Q强度10→0，类别强度保持10。网络/盲物理候选并集、初始类型/Q、Q相关模板形状、固定噪声权重与搜索预算相同。网络Q仍影响候选与模板初始化；类别分数尚未概率校准。

新优化器在[source/pileup/analysis/alpha_gamma_full_prior_fit.py](source/pileup/analysis/alpha_gamma_full_prior_fit.py)。两项强度相等时复现旧算法，旧文件和旧实验不改写。单脉冲搜索参数8 ns、100 ns内近邻使用旧版联合时刻约束，2轮、单个非线性子问题至多20次评价，残差最多补8个分量、加峰阈值9，4 ns候选容量。有限局部搜索不保证全局最优；若干子问题达到评价预算，停止状态原样保存。

## 开发筛选与独立性

[开发协议](protocol.json)在新结果出现前固定：20个已公开开发窗，两γ背景×α源1/10 MHz×各5窗，γ源10 MHz；每窗20 μs，2 ns采样，中央[5,15) μs评价。候选方法只有λQ=0、λC=10，不重新网格调参。

筛选相对物理拟合要求两类总体恢复数不减、逐类完整precision最多下降2个百分点、额外预测≤1.25×原额外数+5；有≥20个目标的条件/粒子组恢复最多下降5个百分点。容差是项目筛选选择，不是文献阈值或统计等效性结论。

{table(dev)}

开发筛选通过：α+34、γ+4、额外预测45→57。γ差的成对窗口95%范围−2.44至+12.25个百分点，尚不稳定。随后[冻结全部设置](frozen-fit-settings.json)，沿用旧U-Net权重及[0.4,0.55]阈值，无重新训练，再生成16个新种子测试窗。旧测试与开发结果不作为新的独立证据。

16窗为两背景×α源1/10 MHz×4窗，γ源10 MHz；种子261091200–261091215，与已有生产种子核对无重用。4个预检窗计入16窗；一个ROOT一个串行G4事件，最多4个独立进程。物理源码、配置和二进制与前轮相同，原生共享SiPM微胞响应不变；完整电压及足够动态范围，不依据ADC削顶拒绝。实际条件见[manifest](manifest.json)，源码另由[native-source-manifest](native-source-manifest.json)核验。

首次启动漏建macros/logs目录，在G4进程启动之前失败，未生成ROOT；补齐后重试，失败日志和说明保存在本地运行记录。所有算法参数、随机种子和物理配置保持冻结。随后外层分析进程以143中断，14窗结果和文件访问记录完整、2窗尚无结果；只对缺少输出的两窗按相同冻结参数与原ROOT补跑，未增加独立样本。首14窗外层返回码未落盘，以结果与访问记录恢复，补跑返回码另记。

## 新独立测试的实际结果

完整恢复：真值Q≥500 p.e.、类型无关一对一匹配≤10 ns、类型正确且|ΔQ/Q|≤20%，中央发射区[5,15) μs。完整precision使用中央重建区的该类预测为分母，分子保留对正电荷保护区/阈下真值的正确类型/Q配对；与目标恢复数分母不同。Q为实现电荷，未标定为沉积能量。

{table(test)}

{decision}

{contrast}

按四源条件分层、2000次成对整窗bootstrap只含有限窗口波动，不含模型/器件系统误差，且每条件仅4窗。γ逐事件相对物理拟合损失{transition['lost']}个原成功事件、补回{transition['gained']}个，不能把净变化当成每事件均改善。弱γ 500–1000 p.e.：{weaktext}。完整分条件、漏检/错分/误拆分/误合并标志、条件Q误差与前史/近邻失败在本地summary及Notebook中，分层有重叠，不作独立因果效应。

## 盲输入、复算与验证

推理接口只接受观测时间/完整电压、冻结标定、网络及其阈值，不接收真实类型、Q、时刻、数量或源率。每个测试窗先完成五种方法，全部16窗结束后才打开标签评分；16个拟合子进程的文件访问检查禁止打开ROOT和truth路径。原始导出和评分由外层程序独立处理。

开发计算在本机临时目录进行，输入来自前轮已归档blind-inputs，并验证哈希和密集图解码。首个真实开发窗的强度0/10端点与旧结果数值一致；重载冻结Torch模型后密集图及端到端五方法也完全一致。3项针对性数值测试通过。实际运行器及验证在`jupyter/run_records/alpha_gamma_full_prior_20261019/`，结果在`jupyter/artifacts/alpha_gamma_full_prior_20261019/`，旧ROOT/模型/结果不覆盖。

需要本地输入归档；独立测试重推理另需`requirements-ml.txt`固定CPU Torch和原U-Net权重。从项目根目录运行，输出路径须不存在：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_full_prior_20261019/replay_analysis.py \\
  --output jupyter/artifacts/alpha_gamma_full_prior_dev_replay --workers 4

bash tools/run_gagg_python.sh pileup/experiments/alpha_gamma_full_prior_20261019/replay_test.py \\
  --output jupyter/artifacts/alpha_gamma_full_prior_test_replay --workers 4
```

`--stems`允许诊断子集，只有完整20/16窗才汇总。便捷包装入口未再整轮重复；本轮实际保存的运行器已完成全部拟合，端点重放、文件访问和重评分范围以本地validation.json为准。重新生成原生样本可从保存的源码在对应Geant4/ROOT环境构建两种二进制，并使用原宏；实际二进制哈希、配置和数据集环境已保存，跨平台重新构建需要另记环境/二进制身份。

## 限制与下一步

同一仿真物理模型下的新种子验证不包括真实器件域偏移。未完成可见电荷到沉积能量的条件标定，条件Q标准差不作为能量分辨率。观测处理耗时分首14窗4并行与补跑2窗2并行记录，每进程1线程，混合包括物理基线与网络；它是离线成本，不作为硬件吞吐率。

继续优先检查弱γ错分、可观测前史及近邻下的类型/Q相关模板与不确定度，在独立标定/开发数据中评价，再决定模板迭代或类型条件软约束。当前另一类Q图未必得到有效监督，不能直接读取两张图充当校准先验。此批新测试今后作为已公开回归，不再用于选参。文献适用条件和审核状态沿用[采用记录L28–L31](../../../docs/literature-use.md#alpha-gamma-full-prior-20261019)，未新增论文性能外推，也未写共享文献库。

过程追加到本地`jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb`第10节，RunInfo/原宏映射在Run_Catalog。本地笔记、模型、ROOT与分析产物按仓库约定不上传。
'''
(folder/'README.md').write_text(body)
anchor='<a id="alpha-gamma-full-prior-20261019"></a>'
readme=project/'README.md';text=readme.read_text();assert anchor not in text
a,g=lookup(test,'class_only','alpha'),lookup(test,'class_only','gamma')
b,h=lookup(test,'baseline','alpha'),lookup(test,'baseline','gamma')
paragraph=f'''{anchor}
2026-10-19 将[仅类别惩罚接回完整混合流程](pileup/experiments/alpha_gamma_full_prior_20261019/README.md)：20个已公开开发窗通过预定逐类恢复/precision/误报筛选后，固定λQ=0、λC=10及旧网络，另生成16个新独立原生窗口。新测试α完整恢复{b['recovered_count']}/{b['target_count']}→{a['recovered_count']}/{a['target_count']}，γ{h['recovered_count']}/{h['target_count']}→{g['recovered_count']}/{g['target_count']}，额外预测{test['diagnostics']['baseline']['unmatched_predictions']}→{test['diagnostics']['class_only']['unmatched_predictions']}。{decision}详见[Notebook第10节](jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)与[L31](docs/literature-use.md#alpha-gamma-full-prior-20261019)。

'''
text=text.replace('## 分析已有数据\n',paragraph+'## 分析已有数据\n',1);readme.write_text(text)
literature=project/'docs/literature-use.md';text=literature.read_text();assert anchor not in text
entry=f'''{anchor}
### L31：2026-10-19——仅类别惩罚的完整混合流程及新独立测试

**依据与边界**：沿用L28–L30的既有精读与原审核状态，没有新增原文复核或论文性能结论。拆分惩罚及筛选容差是项目方案；归一化网络事件分数仍不当作已校准后验。当前Q相关模板仍由初始网络/物理电荷决定，去掉显式Q项不等于去掉全部网络Q影响。

**实现与冻结**：仅λQ 10→0、λC保留10，候选/形状初值、固定噪声权重、时刻/类别/电荷搜索及残差补峰预算沿用旧完整流程。20公开开发窗相对物理拟合α471/791→505/791、γ38/86→42/86，额外预测45→57，达到事先设置的逐类恢复、precision与额外预测筛选。工程容差不作为文献阈值或统计等效性结论。其后冻结旧U-Net权重、阈值[0.4,0.55]、标定和拟合配置，另生成16新种子原生窗，两γ背景×α源1/10 MHz×各4窗，γ源10 MHz。物理与完整电压/足够动态范围前提不变。

**独立计算结果**：同窗物理拟合→仅类别混合，α{b['recovered_count']}/{b['target_count']}→{a['recovered_count']}/{a['target_count']}，γ{h['recovered_count']}/{h['target_count']}→{g['recovered_count']}/{g['target_count']}；额外预测{test['diagnostics']['baseline']['unmatched_predictions']}→{test['diagnostics']['class_only']['unmatched_predictions']}。完整恢复仍要求Q≥500 p.e.、匹配≤10 ns、类型正确、Q误差≤20%，不是沉积能量分辨率。{decision}全部五方法、分条件、错分/漏检与条件Q误差保存在本地summary。{contrast.replace(chr(10), ' ')}

**解释、验证与待办**：γ损失{transition['lost']}个原成功事件、补回{transition['gained']}个；弱γ为{weaktext}。窗口重采样只含有限统计，各分层重叠；若干非线性子问题达到评价预算，不声称全局最优。3项数值检查、旧端点和冻结网络密集图/五方法数值复现、16窗重建文件访问检查、ROOT/RunInfo/原宏与重评分纳入验证。首次缺输出目录在生成ROOT前修复，日志保留；分析外层143中断后仅补跑缺结果的2窗，14个完整结果保持原哈希，不新增样本或调参。下一步在独立标定/开发资料中检验类型/Q/可观测前史条件的模板与不确定度，再设计新实验；此批测试不再用于调参。未验证真实器件、沉积能量与物理域偏移，不修改共享文献库。

**追溯**：[冻结实验](../pileup/experiments/alpha_gamma_full_prior_20261019/README.md)、[Notebook第10节](../jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)、[本地验证](../jupyter/run_records/alpha_gamma_full_prior_20261019/validation.json)。
'''
literature.write_text(text.rstrip()+'\n\n'+entry.rstrip()+'\n')
print('wrote_experiment_readme_and_L31',flush=True)
