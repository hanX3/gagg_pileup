"""Run this frozen intermediate-rate experiment, one ROOT per serial G4 process.

From project root: bash tools/run_gagg_python.sh <this-file> --stage checks
After checking those outputs: <same command> --stage production --workers 2
The original campaign, its ROOTs, and its seeds are preserved.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
spec = importlib.util.spec_from_file_location(
    "frozen_campaign", FOLDER / "source/pileup/analysis/alpha_rate_campaign.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
runner.PROJECT = PROJECT


def save_manifest(manifest):
    temporary = FOLDER / "manifest.json.tmp"
    runner.save(temporary, manifest)
    temporary.replace(FOLDER / "manifest.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("checks", "production"), required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    record = json.loads((FOLDER / "manifest.json").read_text())
    if not 1 <= args.workers <= record["max_parallel_processes"]:
        raise ValueError("Worker count exceeds this experiment's recorded limit")
    if args.stage == "production" and not record["checks_validated"]:
        raise RuntimeError("Validate the small environment checks before production")
    key = "checks" if args.stage == "checks" else "runs"
    plans = record["planned_checks" if key == "checks" else "planned_runs"]
    existing = {r["label"]: r for r in record[key]}
    if any(r["status"] != "completed" for r in existing.values()):
        raise RuntimeError("Inspect incomplete records; they are never silently overwritten")
    pending = [r for r in plans if r["label"] not in existing]
    binary = FOLDER / "build/gagg"
    if runner.digest(binary) != record["binary_sha256"]:
        raise RuntimeError("Binary changed after build validation")
    environment = dict(os.environ)
    environment.update(record["environment"]["g4_dataset_environment"])
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    active = []
    last_started = -float("inf")
    try:
        while pending or active:
            if pending and len(active) < args.workers and time.monotonic() - last_started >= 1.05:
                plan = pending[0]
                output = PROJECT / "g4/data" / (datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss") + ".root")
                if output.exists():
                    time.sleep(1.05)
                    continue
                pending.pop(0)
                macro = FOLDER / "macros" / (plan["label"] + ".mac")
                log = FOLDER / "logs" / (plan["label"] + ".log")
                usage = FOLDER / "logs" / (plan["label"] + ".time.txt")
                if macro.exists() or log.exists() or usage.exists():
                    raise RuntimeError("Run label has existing artifacts")
                macro.write_text(runner.macro(FOLDER, plan["seed"], plan["rate_hz"], output))
                command = ["/usr/bin/time", "-v", "-o", str(usage), str(binary), str(macro)]
                row = dict(plan, sipm_seed=plan["seed"] + 100000, readout_seed=plan["seed"] + 200000,
                           window_us=20., gate_us=[5., 15.], root=str(output.relative_to(PROJECT)),
                           macro=str(macro.relative_to(PROJECT)), macro_sha256=runner.digest(macro),
                           started=datetime.now().astimezone().isoformat(), status="running",
                           command=command, cwd=str((FOLDER / "build").relative_to(PROJECT)),
                           log=str(log.relative_to(PROJECT)), usage_log=str(usage.relative_to(PROJECT)))
                record[key].append(row)
                save_manifest(record)
                stream = log.open("w")
                start = time.monotonic()
                proc = subprocess.Popen(command, cwd=FOLDER / "build", env=environment,
                                        stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                row["pid"] = proc.pid
                save_manifest(record)
                active.append((proc, stream, row, start, output, usage))
                last_started = start
                print("START", row["label"], row["root"], flush=True)
            for task in list(active):
                proc, stream, row, start, output, usage = task
                code = proc.poll()
                if code is None:
                    continue
                stream.close()
                row.update(returncode=code, elapsed_s=time.monotonic() - start,
                           finished=datetime.now().astimezone().isoformat())
                ok = code == 0 and output.exists() and "ROOT file is saved" in (PROJECT / row["log"]).read_text()
                row["status"] = "completed" if ok else "failed"
                if output.exists():
                    row.update(root_bytes=output.stat().st_size, root_sha256=runner.digest(output))
                if usage.exists():
                    match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", usage.read_text())
                    if match:
                        row["max_rss_kb"] = int(match.group(1))
                active.remove(task)
                save_manifest(record)
                print("DONE", row["label"], row["status"], f"{row['elapsed_s']:.1f}s",
                      f"RSS={row.get('max_rss_kb')} KB", flush=True)
                if not ok:
                    raise RuntimeError(f"Failed run; inspect {row['log']}")
            if pending or active:
                time.sleep(.5)
    except BaseException:
        for proc, stream, row, start, output, usage in active:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                proc.wait()
            stream.close()
            row.update(status="interrupted", returncode=proc.returncode, elapsed_s=time.monotonic() - start,
                       finished=datetime.now().astimezone().isoformat())
            if output.exists():
                row["root_bytes"] = output.stat().st_size
        save_manifest(record)
        raise
    print("STAGE COMPLETE", args.stage, len(record[key]), "runs", flush=True)


if __name__ == "__main__":
    main()
