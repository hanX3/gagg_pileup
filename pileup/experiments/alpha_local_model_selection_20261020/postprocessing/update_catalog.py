"""Append actual new ROOT/macro mappings while preserving the old catalog."""
import json
from pathlib import Path

import nbformat


def main():
    project=Path.cwd()
    exp=project/'pileup/experiments/alpha_local_model_selection_20261020'
    record=project/'jupyter/run_records/alpha_local_model_selection_20261020'
    data=json.loads((exp/'manifest.json').read_text())
    assert len(data['runs'])==23 and all(r['status']=='completed' for r in data['runs'])
    path=project/'jupyter/Run_Catalog.ipynb'
    old=nbformat.read(record/'Run_Catalog.ipynb.before',as_version=4)
    current=nbformat.read(path,as_version=4)
    assert current==old,'Catalog changed during production; inspect and retain its new entries first.'
    text='''<a id="alpha-local-model-selection-20261020"></a>
## 2026-10-20：纯α局部脉冲数搜索的独立测试

23个ROOT均在拟合参数冻结后生成，源率0.1/1/10/20/40/60/80/100 MHz，分别8/3/2/2/2/2/2/2窗。每文件一个20 μs WaveformEvent，中央[5,15) μs评价，实际源为p–¹¹B有效单α谱，无γ源。源注入率不等于探测信号率或反应率。完整analog_mV用于盲重建，祖先归属只用于拟合后的评分。

物理/二进制沿用已验证父实验，原生共享微胞历史保留；本次最多3个独立串行G4进程，先验证1/10 MHz两窗再放行其余，预检查包含在23窗中。实际源码、配置、二进制SHA、命令、日志、种子与冻结时序由manifest及freeze记录追溯。ROOT名称只使用秒级时间，不覆盖历史文件；stdout和可执行文件为普通本地临时文件，没有共享盘临时链接。

共享盘mmap读取过慢的末段导出使用普通本机字节副本；生产ROOT仍直接写入并保留在g4/data。使用同一冻结导出器，副本SHA256与原文件一致；60 MHz检验的16个数组逐项相同。21个已启动物理任务均正常结束后恢复剩余只读导出，再启动预定的最后两个100 MHz种子，源码/物理/参数不变；过程保存在manifest.read_cache_transition及本次运行记录。

[实际运行清单](../pileup/experiments/alpha_local_model_selection_20261020/manifest.json) · [参数冻结及科学源码身份](artifacts/alpha_local_model_selection_20261020/freeze.json) · [完整分析结果](artifacts/alpha_local_model_selection_20261020/test-summary.json)

| ROOT | 实际宏 | α注入率/MHz | 源种子 | SiPM种子 | 读出种子 | 状态 | SHA256 |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
'''
    for row in sorted(data['runs'],key=lambda r:r['started']):
        note='（主光子time分支采用冗余字段恢复评分）' if row.get('scoring_field_recovery') else ''
        text+=f"| [{Path(row['root']).name}](../{row['root']}){note} | [{row['label']}](../{row['macro']}) | {row['rate_hz']/1e6:g} | {row['seed']} | {row['sipm_seed']} | {row['readout_seed']} | {row['status']} | `{row['root_sha256']}` |\n"
    text+='''\n一份100 MHz ROOT的PrimaryPhoton时间分支无法解压。评分导出通过唯一的x/y/能量浮点位，与全局WaveformEvent的同光子冗余时间精确一对一连接，覆盖全部光子；另一完整100 MHz文件的16数组等价验证通过，总电压4数组与恢复前相同。原始ROOT及错误记录保留，该分支仍不可直接读取；此适配仅用于评分导出。验证记录为run_records/alpha_local_model_selection_20261020/redundant-time-validation.json。\n'''
    capacity=json.loads((project/data['capacity_extension']['experiment']/'manifest.json').read_text())
    auxiliary=data['failed_attempts']+capacity['validation_runs']
    text+='''\n### 容量上限失败与成对重放检验（不计入23窗分母）

首次100 MHz、20 μs触发max_candidates=10000000保护，控制器随后停止另一进程；两份不完整ROOT保留且不参与拟合。独立容量快照仅将该计算保护提高至30000000，物理/读出概率、C++/二进制和拟合参数不变；相同1 MHz种子的重放验证16个数组逐项一致后重试原100 MHz计划种子。该1 MHz重放只是容量验证，不是额外独立测试。容量扩展及三份附加ROOT的实际配置、日志和状态由对应manifest记录。

| ROOT | 实际宏 | α注入率/MHz | 源种子 | SiPM种子 | 读出种子 | 状态 | SHA256（完整输出） |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
'''
    for row in auxiliary:
        text+=f"| [{Path(row['root']).name}](../{row['root']}) | [{row['label']}](../{row['macro']}) | {row['rate_hz']/1e6:g} | {row['seed']} | {row['sipm_seed']} | {row['readout_seed']} | {row['status']} | `{row.get('root_sha256','不完整，无有效波形')}` |\n"
    cell=nbformat.v4.new_markdown_cell(text,id='alpha-local-production-20261020')
    current.cells.append(cell);nbformat.validate(current);nbformat.write(current,path)
    assert nbformat.read(path,as_version=4).cells[:-1]==old.cells
    (record/'root-macro-index.json').write_text(json.dumps(data['runs'],ensure_ascii=False,indent=2)+'\n')
    (record/'auxiliary-root-macro-index.json').write_text(json.dumps(auxiliary,ensure_ascii=False,indent=2)+'\n')
    print('Catalog: appended 23 test and 3 auxiliary ROOT/macro mappings; old cells preserved.')


if __name__=='__main__':main()
