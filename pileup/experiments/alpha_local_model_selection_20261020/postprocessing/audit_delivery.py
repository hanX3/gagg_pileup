"""Independent checks of the delivered notebooks, artifacts and summary."""
import ast
import base64
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

import nbformat


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    project=Path.cwd()
    record=project/'jupyter/run_records/alpha_local_model_selection_20261020'
    output=project/'jupyter/artifacts/alpha_local_model_selection_20261020'
    exp=project/'pileup/experiments/alpha_local_model_selection_20261020'
    path=project/'jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb'
    old=nbformat.read(record/(path.name+'.before'),as_version=4)
    now=nbformat.read(path,as_version=4)
    nbformat.validate(now)
    old_ids={c.id for c in old.cells}
    assert [c for c in now.cells if c.id in old_ids]==old.cells
    assert len({c.id for c in now.cells})==len(now.cells)
    assert not any(c.id in ('local-alpha-progress-summary','local-alpha-provisional-20261020') for c in now.cells)
    added=[c for c in now.cells if c.id not in old_ids]
    anchors=re.findall(r'<a\s+id="([^"]+)"', '\n'.join(c.source for c in now.cells))
    links=[]
    for c in added:
        if c.cell_type=='code':ast.parse(c.source)
        for href in re.findall(r'\]\(([^)]+)\)', c.source):
            assert not re.search(r'\.(pdf|ipynb|md)(?:#.*)?$',href,re.I),href
            if href.startswith('#'):
                assert anchors.count(href[1:])==1,href
            elif not urlsplit(href).scheme:
                assert (path.parent/unquote(href.split('#')[0])).exists(),href
            links.append(href)
    figures={}
    for number,name in [('17-1','local_controls.png'),('17-2','local_native_rates.png'),('17-3','local_native_examples.png')]:
        assert anchors.count('codex-figure-'+number)==1
        cell=next(c for c in now.cells if c.id=='local-alpha-figure-'+number)
        png=base64.b64decode(cell.outputs[0].data['image/png'])
        assert png==(output/name).read_bytes()
        assert any('图 '+number+'：' in c.source for c in added)
        figures[number]=hashlib.sha256(png).hexdigest()
    catalog_path=project/'jupyter/Run_Catalog.ipynb'
    catalog=nbformat.read(catalog_path,as_version=4)
    old_catalog=nbformat.read(record/'Run_Catalog.ipynb.before',as_version=4)
    nbformat.validate(catalog)
    assert catalog.cells[:-1]==old_catalog.cells
    manifest=json.loads((exp/'manifest.json').read_text())
    summary=json.loads((output/'test-summary.json').read_text())
    assert len(manifest['runs'])==len(summary['windows'])==23
    assert {r['root'] for r in manifest['runs']}=={r['root'] for r in summary['windows']}
    assert {(r['label'],r['rate_hz'],r['seed']) for r in manifest['runs']}=={
        (r['label'],r['rate_hz'],r['seed']) for r in manifest['planned_runs']}
    assert all(r['status']=='completed' for r in manifest['runs'])
    fields=('target_count','recovered_count','unmatched_predictions','missed_count','merge_flags','split_flags')
    for method in ('baseline','local'):
        for rate in summary['rates'][method]:
            rows=[r for r in summary['windows'] if r['rate_hz']==rate['rate_hz']]
            assert len(rows)==rate['windows']
            for key in fields:assert sum(r['methods'][method]['score'][key] for r in rows)==rate[key]
    for row in manifest['runs']:
        assert row['root'] in catalog.cells[-1].source and row['macro'] in catalog.cells[-1].source
        assert row['root_sha256'] in catalog.cells[-1].source
    auxiliary=json.loads((record/'auxiliary-root-macro-index.json').read_text())
    assert len(auxiliary)==3
    for row in auxiliary:
        assert row['root'] in catalog.cells[-1].source and row['macro'] in catalog.cells[-1].source
        assert row['status'] in catalog.cells[-1].source
    assert json.loads((record/'capacity-equivalence.json').read_text())['passed']
    validation=json.loads((record/'validation.json').read_text())
    assert validation['passed'] and validation['scoring_comparisons']==46
    document=project/'docs/literature-use.md'
    assert document.read_bytes().startswith((record/'literature-use.md.before').read_bytes())
    assert document.read_text().count('### L40：2026-10-20')==1
    frozen=json.loads((output/'freeze.json').read_text())
    for name,expected in frozen['frozen_files'].items():assert sha(project/name)==expected,name
    source=json.loads((record/'postprocessing-source.json').read_text())
    for name,expected in source.items():assert sha(project/name)==expected,name
    result=dict(passed=True,original_notebook_cells=len(old.cells),notebook_cells=len(now.cells),
        added_cells=len(added),all_old_cells_outputs_exactly_equal=True,
        catalog_old_cells_exactly_equal=True,new_root_macro_mappings=23,
        auxiliary_root_macro_mappings=3,capacity_extension_paired_validation_passed=True,
        planned_runs_exactly_matched=True,new_figure_sha256=figures,
        new_relative_and_internal_links_checked=len(links),
        all_eight_rate_summaries_independently_reaggregated=True,
        frozen_source_files_unchanged=len(frozen['frozen_files']),
        postprocessing_source_files_checked=len(source),
        notebook_sha256=sha(path),catalog_sha256=sha(catalog_path),
        test_summary_sha256=sha(output/'test-summary.json'))
    (record/'delivery-audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
