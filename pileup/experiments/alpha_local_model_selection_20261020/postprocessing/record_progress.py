"""Make the completed test subset visible while high-rate production runs."""
from datetime import datetime
import json
from pathlib import Path
import nbformat

p=Path.cwd();out=p/'jupyter/artifacts/alpha_local_model_selection_20261020';rec=p/'jupyter/run_records/alpha_local_model_selection_20261020'
manifest=json.loads((p/'pileup/experiments/alpha_local_model_selection_20261020/manifest.json').read_text())
rows=[json.loads(f.read_text()) for f in (out/'test/selected').glob('*.json')]
lines=['| α源注入率/MHz | 已评分窗口 | 目标数 | 原基线完整恢复 | 局部搜索完整恢复 | 原→新额外预测 |',
       '| ---: | ---: | ---: | --- | --- | --- |']
for rate in sorted({r['rate_hz'] for r in rows}):
    rr=[r for r in rows if r['rate_hz']==rate]
    n=sum(r['methods']['baseline']['score']['target_count'] for r in rr)
    counts={m:sum(r['methods'][m]['score']['recovered_count'] for r in rr) for m in ('baseline','local')}
    extra={m:sum(r['methods'][m]['score']['unmatched_predictions'] for r in rr) for m in ('baseline','local')}
    def fmt(m):return f"{counts[m]}/{n}（{100*counts[m]/n:.1f}%）" if n else '0/0'
    lines.append(f"| {rate/1e6:g} | {len(rr)} | {n} | {fmt('baseline')} | {fmt('local')} | {extra['baseline']}→{extra['local']} |")
table='\n'.join(lines);now=datetime.now().astimezone().isoformat()
save=dict(recorded_at=now,scored_windows=len(rows),planned_windows=23,rows=rows,table_markdown=table,
          native_completed=sum(r['status']=='completed' for r in manifest['runs']))
(out/'interim-summary.json').write_text(json.dumps(save,ensure_ascii=False,indent=2)+'\n')
path=p/'jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb'
before=nbformat.read(rec/'codex_0003_Alpha_Pileup_Rate_Fitting.ipynb.before',as_version=4)
current=nbformat.read(path,as_version=4)
provisional=[c for c in current.cells if c.id in ('local-alpha-progress-summary','local-alpha-provisional-20261020')]
if provisional:
    previous=rec/'notebook-provisional-17-scored.ipynb'
    if not previous.exists():nbformat.write(current,previous)
    current.cells=[c for c in current.cells if c.id not in ('local-alpha-progress-summary','local-alpha-provisional-20261020')]
assert current==before
notice=nbformat.v4.new_markdown_cell('''**2026-10-20计算进度**：[第17节](#pure-alpha-local-model-selection-20261020)已实现并列出当前完成的新原生同窗比较。高率生产仍进行，表格明确标为阶段结果；完整比较和图17-1至17-3会在同一章节更新。原第14节36窗为历史基线。\n''',id='local-alpha-progress-summary')
text=f'''<a id="pure-alpha-local-model-selection-20261020"></a>
## 17. 主动拆峰、局部脉冲数搜索与多初值联合拟合（计算进行中）

**进度记录：{now}。当前仅{len(rows)}/23个新原生窗口完成两方法评分，不是完整八档扫描结果。** 原生已完成{save['native_completed']}/23窗；剩余高率生产和导出继续运行。最终图号预留17-1、17-2、17-3，尚未在本节生成。

新实现保留95/351 ns固定双衰减模板，盲估计K/t/Q；增加主动单→双拆分、残差补峰、删除/合并及多初值局部联合调时，超过64分量也允许局部联合更新。独立48窗单α标定、24窗开发单α及三个公开1/10/40 MHz窗选择参数，17:34:07 +08:00冻结后才生成新测试。拟合不接收真实K、时间、电荷和源率；真值只在推断后评分。完整动态范围、原生共享微胞、2 ns采样及0.1 mV电子噪声保持。

### 17.1 当前已完成的新测试同窗比较

{table}

完整恢复：真实Q≥500 p.e.、一对一时间匹配±10 ns、电荷误差≤20%。电荷不是已标定沉积能量。阶段结果保留低率改善与退步，不据小样本定义最大计数率或分辨极限。

### 17.2 已完成的检查与后续

五项数值检查通过，覆盖梯度、空/单脉冲、隐藏双脉冲、冷/热启动一致性及超过64分量的局部更新。174条同均值模板软件控制（仅电子白噪声）中，空/单无额外预测；126条双脉冲的两个事件均恢复且无额外输出案例84→94。该控制不是实物时间分辨率。全部23窗完成后，统一补齐逐率恢复、额外预测、漏检、误拆分/合并、条件电荷/时间误差及时间错位对照，不用新测试重新选参。

[阶段结果](artifacts/alpha_local_model_selection_20261020/interim-summary.json) · [冻结参数与源码身份](artifacts/alpha_local_model_selection_20261020/freeze.json) · [实际原生运行清单](../pileup/experiments/alpha_local_model_selection_20261020/manifest.json)
'''
current.cells.insert(1,notice);current.cells.append(nbformat.v4.new_markdown_cell(text,id='local-alpha-provisional-20261020'))
nbformat.validate(current);nbformat.write(current,path);nbformat.write(current,rec/'notebook-provisional.ipynb')
old_ids={c.id for c in before.cells};assert [c for c in current.cells if c.id in old_ids]==before.cells
print('Provisional Notebook section recorded with',len(rows),'scored windows; all old cells preserved.')
