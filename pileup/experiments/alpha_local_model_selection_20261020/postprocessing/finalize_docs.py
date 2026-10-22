"""Append concise result pointers without replacing existing README content."""
import hashlib
import json
from pathlib import Path


def main():
    p=Path.cwd();r=p/'jupyter/run_records/alpha_local_model_selection_20261020'
    out=p/'jupyter/artifacts/alpha_local_model_selection_20261020'
    exp=p/'pileup/experiments/alpha_local_model_selection_20261020'
    summary=json.loads((out/'test-summary.json').read_text())
    report=json.loads((out/'report.json').read_text())
    frozen=json.loads((out/'freeze.json').read_text())
    manifest=json.loads((exp/'manifest.json').read_text())
    assert len(summary['windows'])==23 and all(x['status']=='completed' for x in manifest['runs'])
    index={m:{x['rate_hz']:x for x in summary['rates'][m]} for m in ('baseline','local')}
    def counts(rate):
        a,b=(index[m][rate] for m in ('baseline','local'))
        return f"{a['recovered_count']}/{a['target_count']}→{b['recovered_count']}/{b['target_count']}"
    def extras(rate):return '→'.join(str(index[m][rate]['unmatched_predictions']) for m in ('baseline','local'))
    root=f'''

2026-10-20 实现[纯α主动拆峰和局部模型搜索](pileup/experiments/alpha_local_model_selection_20261020/README.md)：固定95/351 ns模板，未知K/t/Q，增加单→双拆分、增删/合并及多初值联合重拟合。独立标定/开发后冻结，在23个新原生窗口上与原基线同窗比较，20/40 MHz完整恢复分别{counts(2e7)}、{counts(4e7)}，10 MHz为{counts(1e7)}；40 MHz额外预测{extras(4e7)}。八档完整表、误差、失败和时间错位对照见[本地Notebook0003第17节](jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-local-model-selection-20261020)及[采用记录L40](docs/literature-use.md#pure-alpha-local-model-selection-20261020)。低残差不能替代逐事件恢复；高率每档仅2窗，电荷尚未连接沉积能标。
'''
    pileup=f'''

## 纯α的主动拆峰与局部模型搜索（2026-10-20）

[实现模块](analysis/alpha_local_fit.py)在原盲拟合初值上增加单→双拆分、残差补峰、删除/合并和多初值局部联合时间/非负幅度优化，超过64分量也执行局部联合更新。固定95/351 ns和慢权重，完整analog_mV输入；未知真实K/t/Q，不输入源率。48个独立单α标定窗、24个开发单α及三个公开率窗选定经验SSE代价后，冻结并生成23个新原生测试窗。

0.1–100 MHz八档同窗比较保留改善和退步：20/40 MHz完整恢复{counts(2e7)}、{counts(4e7)}，10 MHz为{counts(1e7)}；40 MHz额外预测{extras(4e7)}。174条固定均值＋白噪声控制单列；高率两窗不能确定最高可用率，Q恢复不等于沉积能量恢复。新增[冻结快照与重跑入口](experiments/alpha_local_model_selection_20261020/README.md)、[本地结果与固定图17-1至17-3](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-local-model-selection-20261020)和[L40来源/边界](../docs/literature-use.md#pure-alpha-local-model-selection-20261020)，原36窗基线和旧输出保留。
'''
    experiment=f'''

## 本轮完成结果

参数于{frozen['frozen_at']}冻结；选择倍率{frozen['selected']['factor']:g}，电荷项系数为{frozen['selected']['config']['addition_cost_mV2_per_pe']:.8f} mV²/p.e.。23份原生ROOT及两方法拟合全部完成；20/40/100 MHz完整恢复分别为{counts(2e7)}、{counts(4e7)}、{counts(1e8)}，10 MHz为{counts(1e7)}。40/100 MHz额外预测为{extras(4e7)}、{extras(1e8)}；不能只按净恢复数评价。匹配后的电荷误差、同一目标新增/丢失、漏检与拆分/合并标志、错位时间对照保存在本地test-summary.json和report.json。

只读导出的I/O调整见manifest.read_cache_transition：生产ROOT始终直接写在共享g4/data；经60 MHz文件16个数组逐项等价检查后，末段使用校验过SHA256的普通/tmp字节副本读取。原始ROOT、未完成导出的中间文件和日志均保留，拟合/物理/冻结参数没有改变。该调度与I/O改动的实际代码保存在postprocessing/，推断仍使用inference-source/。

五项数值检查、无真值输入冷重放、23个原生ROOT/宏/RunInfo/种子/冻结时序、46组模型重构和重评分、旧Notebook输出保留由本地run_records/alpha_local_model_selection_20261020的验证JSON追溯。绘图/交付脚本为postprocessing/alpha_local_fit_report.py及对应可复用模块；从项目根目录仅重画三张图使用：

```bash
env MPLCONFIGDIR=/tmp/gagg-alpha-local-report \\
  bash tools/run_gagg_python.sh -m pileup.analysis.alpha_local_fit_report --figures-only
```

完整的文档交付入口对旧笔记/文档有保护检查，用于本轮一次性交付；已交付后使用--figures-only或--validate-only。重拟合使用前述新输出目录命令。高率每档两窗、局部搜索预算、固定平均形状失配及偶然匹配限制均保留；不据此给出普适时间分辨率、最大计数率、实物或实时吞吐结论。

100 MHz首次触发max_candidates=10000000计算保护，原失败ROOT及日志保留在manifest.failed_attempts；独立alpha_local_capacity_20261020配置快照仅将该保护上限提高至30000000。相同1 MHz种子的16数组逐项一致验证通过后，重试原计划两套100 MHz种子；拟合、物理/读出概率和二进制保持。该配对重放不计入23窗测试分母，三份附加ROOT/宏映射也登记在本地Run_Catalog。

一份100 MHz ROOT的主光子到达时间分支无法解压，原文件保留且该分支仍不可直接读取。评分导出用全局WaveformEvent中的冗余时间，按唯一的x/y/能量浮点位精确一对一连接并核对完整覆盖；另一完整100 MHz文件16个导出数组逐项一致，恢复前后总波形4数组也相同。恢复代码为postprocessing/root_primary_time_recovery.py及recover_primary_times.py，验证在本地redundant-time-validation.json；冻结拟合源码/参数和总电压输入不变。重放拟合直接使用已验证的缓存；重新导出该文件评分字段需启用此适配器。
'''
    changed=[]
    for path,addition,backup_name in (
        (p/'README.md',root,'root_README.md.before_final_report'),
        (p/'pileup/README.md',pileup,'pileup_README.md.before_final_report'),
        (exp/'README.md',experiment,'experiment_README.md.before_final_report')):
        original=path.read_bytes();backup=r/backup_name;assert not backup.exists()
        assert addition.strip() not in path.read_text()
        backup.write_bytes(original);path.write_bytes(original+addition.encode())
        assert path.read_bytes().startswith(original)
        changed.append(dict(path=str(path.relative_to(p)),old_prefix_preserved=True,
            old_sha256=hashlib.sha256(original).hexdigest(),new_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (r/'readme-delivery.json').write_text(json.dumps(dict(passed=True,changed=changed),ensure_ascii=False,indent=2)+'\n')
    print('Appended final result pointers to three READMEs; prior contents preserved.')


if __name__=='__main__':main()
