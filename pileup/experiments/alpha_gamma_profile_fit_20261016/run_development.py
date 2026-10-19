from pathlib import Path
import sys,runpy
folder=Path(__file__).resolve().parent
phase=int(sys.argv.pop(1));assert phase in (1,2,3,4)
module=sys.argv.pop(1);assert module.startswith('pileup.analysis.')
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
override=folder/f'development-source-v{phase}'
if phase==4:override=override/'pileup/analysis'
pileup.analysis.__path__.insert(0,str(override))
runpy.run_module(module,run_name='__main__')
