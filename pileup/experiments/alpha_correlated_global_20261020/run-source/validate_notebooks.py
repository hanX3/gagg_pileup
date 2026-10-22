"""Validate only the new notebook work and preservation of prior content."""
import base64
import json
from pathlib import Path
import re

import nbformat

project=Path(__file__).resolve().parents[3]
record=project/'jupyter/run_records/alpha_correlated_global_20261020'
artifacts=project/'jupyter/artifacts/alpha_correlated_global_20261020'
notebook=nbformat.read(project/'jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb',as_version=4)
original=nbformat.read(record/'before__jupyter__codex_0003_Alpha_Pileup_Rate_Fitting.ipynb.before',as_version=4)
catalog=nbformat.read(project/'jupyter/Run_Catalog.ipynb',as_version=4)
old_catalog=nbformat.read(record/'before__jupyter__Run_Catalog.ipynb.before',as_version=4)
nbformat.validate(notebook);nbformat.validate(catalog)
assert notebook.cells[:94]==original.cells and len(original.cells)==94
assert catalog.cells[:len(old_catalog.cells)]==old_catalog.cells
new=notebook.cells[94:]
markdown='\n'.join(c.source for c in new if c.cell_type=='markdown')
for number in range(1,6):
    assert markdown.count(f'id="codex-figure-19-{number}"')==1
    assert f'0003-19-{number}' in markdown
for cells in (new,catalog.cells[len(old_catalog.cells):]):
    for cell in cells:
        if cell.cell_type=='markdown':
            links=re.findall(r'\[[^\]]*\]\(([^)]+)\)',cell.source)
            assert all(not re.search(r'\.(?:pdf|ipynb|md)(?:#.*)?$',link,re.I) for link in links),links
        else:
            assert cell.execution_count is not None and all(o.output_type!='error' for o in cell.outputs)
images=[o.data['image/png'] for c in new if c.cell_type=='code' for o in c.outputs if o.output_type=='display_data' and 'image/png' in o.data]
assert len(images)==5,len(images)
for encoded,filename in zip(images,['covariance_calibration.png','native_recovery.png','native_diagnostics.png','software_controls.png','native_examples.png']):
    assert base64.b64decode(encoded)==(artifacts/filename).read_bytes(),filename
result=dict(nbformat_passed=True,old_0003_cells_unchanged=94,
    old_0003_code_cells_preserved=sum(c.cell_type=='code' for c in original.cells),
    old_0003_outputs_preserved=sum(len(c.get('outputs',[])) for c in original.cells),
    catalog_previous_cells_unchanged=len(old_catalog.cells),new_figures_embedded_and_equal_to_artifacts=5,
    figure_ids=[f'0003-19-{i}' for i in range(1,6)],new_forbidden_links=0,
    appended_0003_code_cells_executed=sum(c.cell_type=='code' for c in new),
    historical_cells_rerun=False)
(record/'notebook-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print('NOTEBOOK_VALIDATION',result)
