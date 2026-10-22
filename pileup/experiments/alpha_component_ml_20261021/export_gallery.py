"""Replay archived gallery math with a fixed legend position for rendering.

An earlier long-lived coordinator terminated with exit 143 while exporting
figures. Partial outputs are archived, not deleted. Frozen predictions, case
selection, component synthesis, all samples and figure IDs remain unchanged.
"""
import json
import os
from pathlib import Path
import sys

PROJECT=Path(__file__).resolve().parents[3]
EXP=Path(__file__).resolve().parent
OUT=PROJECT/'jupyter/artifacts/alpha_component_ml_20261021'
REC=PROJECT/'jupyter/run_records/alpha_component_ml_20261021'
os.environ['GAGG_COMPONENT_PROJECT']=str(PROJECT)
sys.path.insert(0,str(EXP/'inference-source'))

from pileup.analysis.alpha_component_evaluate import check_freeze
from pileup.analysis.alpha_rate_campaign import digest,save
from pileup.analysis import alpha_component_report as report
import matplotlib


def main():
    check_freeze()
    old=OUT/'gallery';archive=REC/'partial-gallery-pass1'
    assert old.exists() and not archive.exists()
    preserved=[dict(path=str(p.relative_to(old)),bytes=p.stat().st_size,sha256=digest(p))
               for p in old.rglob('*') if p.is_file()]
    old.rename(archive)
    save(REC/'gallery-rendering-replay.json',dict(previous_exit_code=143,
         cause='No Python traceback; cause unconfirmed. Coordinator terminated during gallery export.',
         retained_partial_outputs=str(archive.relative_to(PROJECT)),retained_files=preserved,
         frozen_source_changed=False,predictions_changed=False,figure_ids_changed=False,
         rendering_change='matplotlib legend.loc fixed to upper right, avoids automatic best-position scan',
         adapter_sha256=digest(__file__)))
    matplotlib.rcParams['legend.loc']='upper right'
    report.gallery()
    print('GALLERY_EXPORT_COMPLETE',flush=True)


if __name__=='__main__':main()
