from pathlib import Path
import json
import runpy
import sys
access=[]
def audit(event,args):
    if event=='open' and isinstance(args[0],(str,bytes)):
        path=str(args[0])
        if path.endswith('.root') or '/truth/' in path:
            raise RuntimeError('Blind fit attempted ROOT/truth access: '+path)
        access.append(path)
sys.addaudithook(audit)
stem=sys.argv[sys.argv.index('--stem')+1]
runpy.run_path(str(Path(__file__).with_name('run_test.py')),run_name='__main__')
path=Path(__file__).parent/'test-logs'/(stem+'.file-access.json')
path.write_text(json.dumps(dict(stem=stem,root_truth_guard_enabled=True,fit_completed=True,
    opened_files=sorted(set(access)),number_of_opens=len(access)),indent=2)+'\n')
