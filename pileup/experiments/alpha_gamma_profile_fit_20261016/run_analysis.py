from pathlib import Path
import sys,runpy
folder=Path(__file__).resolve().parent
sys.path.insert(0,str(folder/'source'))
module=sys.argv.pop(1)
assert module.startswith('pileup.analysis.')
runpy.run_module(module,run_name='__main__')
