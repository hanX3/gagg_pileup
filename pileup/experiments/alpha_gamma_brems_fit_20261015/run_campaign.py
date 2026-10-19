from pathlib import Path
import importlib.util
import sys
folder = Path(__file__).resolve().parent
project = folder.parents[2]
sys.path.insert(0, str(project))
sys.path.insert(0, str(folder/"source"))
spec = importlib.util.spec_from_file_location("frozen_mixed_campaign", folder/"source/pileup/analysis/alpha_gamma_campaign.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.PROJECT = project
module.main(default_folder=folder)
