from pathlib import Path
import hashlib
import json
import sys
import nbformat as nb

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
work=Path('/tmp/gagg-full-prior-20261019')
sys.path.insert(0,str(work/'source'));sys.path.insert(0,str(work/'dependency-source'))
from pileup.analysis.run_catalog import collect
relative=Path('pileup/experiments/alpha_gamma_full_prior_20261019/manifest.json')
manifest=json.loads((project/relative).read_text())
assert len(manifest['runs'])==16 and all(r['status']=='completed' for r in manifest['runs'])
rows=collect(project,only_roots=[r['root'] for r in manifest['runs']],campaign_manifests=[relative])
assert len(rows)==16
for row in rows:
    assert row['waveform_events']==1 and row['window_us']==20.
    row['role']='test_mixed'
(record/'run-catalog-rows.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
path=project/'jupyter/Run_Catalog.ipynb'
book=nb.read(path,as_version=4);before=nb.read(record/'Run_Catalog.ipynb.before',as_version=4)
assert book.cells==before.cells
anchor='alpha-gamma-full-prior-native-20261019'
assert not any(anchor in c.get('source','') for c in book.cells)
cells=[nb.v4.new_markdown_cell('''<a id="alpha-gamma-full-prior-native-20261019"></a>
## 2026-10-19：α＋γ完整混合流程电荷/类别惩罚的独立测试

16个新种子原生窗口：两种γ背景×α源1/10 MHz×各4窗，γ源10 MHz；20 μs完整窗、中央[5,15) μs评价。
4个预检窗计入16窗，一个ROOT仅一个串行G4事件，同时最多4个独立任务；宏与ROOT核对。
物理源码、SiPM配置和二进制与前轮相同，采用完整电压及足够动态范围，不依据ADC削顶拒绝。
新种子261091200–261091215，算法/模型/阈值先冻结再生成，旧ROOT和旧结果保留。
首次启动漏建宏/日志目录时尚未启动G4，未产生ROOT；补齐后重试，初次失败日志保留。

[实际manifest](../pileup/experiments/alpha_gamma_full_prior_20261019/manifest.json)、
[冻结设置](../pileup/experiments/alpha_gamma_full_prior_20261019/frozen-fit-settings.json)、
[分析Notebook第10节](codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)。
以下条件从实际RunInfo读取；原宏及种子映射保留。
'''),nb.v4.new_code_cell('''from pathlib import Path
import json
import pandas as pd
_full_catalog_project=next(p for p in [Path.cwd(),*Path.cwd().parents] if (p/'g4/data').is_dir())
_full_catalog_rows=json.loads((_full_catalog_project/'jupyter/run_records/alpha_gamma_full_prior_20261019/run-catalog-rows.json').read_text())
pd.DataFrame(_full_catalog_rows)[['root_file','role','status','gamma_model','alpha_injection_rate_Hz',
    'gamma_injection_rate_Hz','seed','waveform_events','window_us','macro_file','notes']]
''')]
import ast
namespace={};tree=ast.parse(cells[1].source);last=tree.body.pop()
exec(compile(tree,'<full-prior-catalog>','exec'),namespace)
value=eval(compile(ast.Expression(last.value),'<full-prior-catalog>','eval'),namespace)
count=max((c.get('execution_count') or 0 for c in before.cells if c.cell_type=='code'),default=0)+1
cells[1].execution_count=count
cells[1].outputs=[nb.v4.new_output('execute_result',data={'text/plain':repr(value),'text/html':value._repr_html_()},execution_count=count)]
book.cells.extend(cells);nb.validate(book);nb.write(book,path)
assert nb.read(path,as_version=4).cells[:-2]==before.cells
(record/'run-catalog-validation.json').write_text(json.dumps(dict(roots=16,original_cells=len(before.cells),
    appended_cells=2,old_cells_preserved=True,runinfo_checked=True,
    notebook_sha256=hashlib.sha256(path.read_bytes()).hexdigest()),indent=2)+'\n')
print('catalog_appended',len(rows),'original_cells',len(before.cells),flush=True)
