from pathlib import Path
import hashlib
import json
import shutil
import nbformat

project=Path(__file__).resolve().parents[3];record=Path(__file__).resolve().parent
path=project/'jupyter/Run_Catalog.ipynb'
book=nbformat.read(path,as_version=4)
old=nbformat.read(record/'Run_Catalog.ipynb.before',as_version=4)
assert book.cells[:-2]==old.cells
note='分析阶段外层进程以143中断，已保留14窗全部盲结果；仅对无结果的后2窗按相同冻结配置和原ROOT补跑，不增加独立样本。首14窗外层返回码未落盘，以完整结果与文件访问记录恢复；2窗补跑返回码另记。该分析中断不改变16个ROOT的原生有效状态。'
assert note not in book.cells[-2].source
outputs=book.cells[-1].outputs
book.cells[-2].source+='\n\n'+note+'\n'
nbformat.validate(book);nbformat.write(book,path)
assert book.cells[-1].outputs==outputs and book.cells[:-2]==old.cells
shutil.copy2(record/'run-catalog-validation.json',record/'run-catalog-validation-before-analysis-interruption.json')
value=json.loads((record/'run-catalog-validation.json').read_text())
value.update(analysis_interruption_not_a_native_data_failure=True,notebook_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(record/'run-catalog-validation.json').write_text(json.dumps(value,indent=2)+'\n')
print('catalog_analysis_interruption_noted_old_cells_outputs_preserved',flush=True)
