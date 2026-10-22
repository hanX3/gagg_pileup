from pathlib import Path
import sys
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'inference-source'))
from pileup.analysis import alpha_local_fit_study as module
module.PROJECT=folder.parents[2]
module.EXPERIMENT=folder
module.OUTPUT=module.PROJECT/'jupyter/artifacts/alpha_local_model_selection_20261020'
module.RECORD=module.PROJECT/'jupyter/run_records/alpha_local_model_selection_20261020'
if __name__=='__main__':module.main()
