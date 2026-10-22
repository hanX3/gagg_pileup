"""Four independent fitting processes after native production, on this 4-CPU VM.

Scheduling only. Frozen inference and all completed window/method caches are
retained. Earlier two schedulers were interrupted to reduce resource contention.
"""
import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess
import sys
import time

project = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project))
from pileup.analysis.alpha_correlated_evaluate import frozen, physical_window, summarize_test
from pileup.analysis.alpha_correlated_study import EXPERIMENT, OUTPUT, RECORD
from pileup.analysis.alpha_rate_campaign import save


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--root'); a = parser.parse_args()
    if a.root:
        config = json.loads((OUTPUT/'freeze.json').read_text())
        manifest = json.loads((EXPERIMENT/'manifest.json').read_text())
        row = next(r for r in manifest['runs'] if Path(r['root']).stem == a.root)
        assert row['status'] == 'completed' and 'export_metadata' in row
        physical_window(row, config['selected'])
        return
    print('WAIT_FOR_NATIVE_COMPLETION; inference unchanged', flush=True)
    last_native = -1
    while True:
        manifest = json.loads((EXPERIMENT/'manifest.json').read_text())
        assert all(r['status'] not in ('failed', 'interrupted') for r in manifest['runs'])
        complete = sum(r['status'] == 'completed' and 'export_metadata' in r for r in manifest['runs'])
        if complete != last_native:
            print('NATIVE_EXPORTS_READY', complete, '/23', flush=True); last_native = complete
        if complete == 23:
            break
        time.sleep(2.)
    frozen()
    active = {}; history = []
    rows = manifest['runs']; finished = {}
    for row in rows:
        path = OUTPUT/'test/windows'/(Path(row['root']).stem+'.json')
        if path.exists(): finished[row['root']] = json.loads(path.read_text())
    pending = [r for r in rows if r['root'] not in finished]
    try:
        while pending or active:
            while pending and len(active) < 4:
                row = pending.pop(0); stem = Path(row['root']).stem
                command = [sys.executable, str(Path(__file__).resolve()), '--root', stem]
                start = datetime.now().astimezone().isoformat()
                process = subprocess.Popen(command, cwd=project)
                active[stem] = (process, command, start, time.monotonic(), row)
                print('AFTER_NATIVE_FIT_START', stem, flush=True)
            for stem, (process, command, start, tick, row) in list(active.items()):
                status = process.poll()
                if status is None: continue
                history.append(dict(root_stem=stem, command=command, started=start,
                                    elapsed_s=time.monotonic()-tick, returncode=status))
                save(RECORD/'after-native-process-evaluation.json', dict(workers=4, inference_changed=False,
                    earlier_completed_outputs_reused=True, starts_after_native_completion=True, jobs=history))
                assert status == 0, (stem, status)
                del active[stem]
                path = OUTPUT/'test/windows'/(stem+'.json')
                finished[row['root']] = json.loads(path.read_text())
                print('AFTER_NATIVE_FIT_DONE', len(finished), '/23', flush=True)
            time.sleep(.5)
    finally:
        for process, command, start, tick, row in active.values():
            if process.poll() is None:
                process.terminate(); process.wait()
    assert len(finished) == 23
    summarize_test(list(finished.values()))
    print('ALL_PHYSICAL_DONE', flush=True)


if __name__ == '__main__': main()
