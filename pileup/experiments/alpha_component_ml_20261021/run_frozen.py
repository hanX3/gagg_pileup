"""Execute the initial production/analysis using only saved Python sources."""
from pathlib import Path
import argparse
import os
import sys

project=Path(__file__).resolve().parents[3]
os.environ['GAGG_COMPONENT_PROJECT']=str(project)
sys.path.insert(0,str(Path(__file__).resolve().parent/'inference-source'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('native','evaluate','controls','diagnostics','figures','gallery','validate'))
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    from pileup.analysis.alpha_component_evaluate import check_freeze
    check_freeze()
    if a.stage=='native':
        from pileup.analysis.alpha_component_study import native
        native(a.workers)
    elif a.stage in ('evaluate','controls','diagnostics'):
        from pileup.analysis import alpha_component_evaluate as evaluation
        evaluation.evaluate(a.workers) if a.stage=='evaluate' else getattr(evaluation,a.stage)()
    elif a.stage in ('figures','gallery'):
        from pileup.analysis import alpha_component_report as report
        getattr(report,a.stage)()
    else:
        from pileup.analysis.alpha_component_validate import validate
        validate()


if __name__=='__main__':main()
