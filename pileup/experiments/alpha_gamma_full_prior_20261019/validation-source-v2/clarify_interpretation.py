"""Add the observed classification tradeoff without changing code or outputs."""
from pathlib import Path
import json
import nbformat

project=Path(__file__).resolve().parents[3]
record=Path(__file__).resolve().parent
path=project/'jupyter/codex_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb'
book=nbformat.read(path,as_version=4);before=nbformat.from_dict(json.loads(nbformat.writes(book)))
note='''补充分型诊断：相对物理拟合损失的8个原成功γ全部被错分为α，同时补回7个，净少1个。弱γ 500–1000 p.e.完整恢复16/52→13/52，错类型22→29。同批完整混合从无显式惩罚(0,0)到仅类别(0,10)，α367/643→382/643、γ32/79→29/79、额外预测49→54；其中5个原成功γ转为错类型、补回2个。去掉Q项相对旧两项混合则γ28/79→29/79。类别提示也有取舍，不能把全部失败归给显式Q项。

相对物理拟合的γ恢复率差为−1.27个百分点，成对整窗95%范围[−6.06,+3.70]个百分点，包含零，不宣称总体显著下降；高α率轫致辐射组5/28→1/28和低α率该背景的α24/29→21/29则触发了预定筛选失败。下一步在训练/独立标定资料中优先做Q匹配的独立PSD消融、类别分数校准，以及可观测前史/Q条件模板检查，不用本批测试继续选参。'''
target=next(c for c in book.cells[49:] if c.cell_type=='markdown' and c.source.startswith('### 10.3 '))
assert '补充分型诊断' not in target.source;target.source+='\n\n'+note
nbformat.validate(book);nbformat.write(book,path)
after=nbformat.read(path,as_version=4)
assert after.cells[:49]==before.cells[:49]
assert [c for c in after.cells if c.cell_type=='code']==[c for c in before.cells if c.cell_type=='code']
short='同批完整混合由无显式惩罚到仅类别提示，α367/643→382/643、γ32/79→29/79，额外预测49→54；丢失5个原成功γ且均错分α，补回2个。相对物理拟合损失的8个原成功γ也全部错分α。类别提示仍有取舍，下一步优先在训练/独立标定资料中做Q匹配PSD、类别分数校准和可观测前史/Q条件模板检查，不用本批测试选参。γ相对物理拟合的窗口区间跨零，不声称总体显著下降。'
p=project/'pileup/experiments/alpha_gamma_full_prior_20261019/README.md'
text=p.read_text();text=text.replace('## 限制与下一步\n','## 补充分型诊断\n\n'+short+'\n\n## 限制与下一步\n',1);p.write_text(text)
p=project/'docs/literature-use.md';p.write_text(p.read_text().rstrip()+'\n\n**补充分型对照**：'+short+'\n')
(record/'interpretation-update.json').write_text(json.dumps(dict(original_49_cells_preserved=True,
    all_code_cells_and_outputs_preserved=True,notebook_format_valid=True,algorithm_or_metric_changes=False),indent=2)+'\n')
print('interpretation_clarified_all_code_outputs_preserved',flush=True)
