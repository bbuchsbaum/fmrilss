#!/usr/bin/env python3
"""Rerun the review checks serially; preserve the recorded_results directory."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
scripts = [
    ['tmp/glmsingle/check_algebra.py'],
    ['tmp/glmsingle/review_precision/run_reviewer_checks.py'],
    ['tmp/glmsingle/review_precision/check_precision.py'],
    ['tmp/glmsingle/source_edge_probe.py', '--output', 'tmp/glmsingle/due_diligence/source_edge_results.json'],
    ['tmp/glmsingle/due_diligence/pc_prefix_benchmark.py', '--output', 'tmp/glmsingle/due_diligence/pc_prefix_results.json'],
]
logs = root/'rerun_logs'
logs.mkdir(exist_ok=True)
for args in scripts:
    print('Running', args[0], flush=True)
    proc = subprocess.run([sys.executable, *args], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (logs/(Path(args[0]).stem+'.txt')).write_text(proc.stdout)
    if proc.returncode:
        print(proc.stdout)
        raise SystemExit(proc.returncode)
    print('  completed; log:', str((logs/(Path(args[0]).stem+'.txt')).relative_to(root)), flush=True)
print('All checks completed. Recorded results are unchanged.')
