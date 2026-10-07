#!/usr/bin/env bash
# Run the fmrilss vs Python LSS comparison end to end.
# Usage: bash run_all.sh <work_dir>
# Requires R with fmrilss installed, and Python with numpy and nilearn.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
work="${1:-sim_runs}"
mkdir -p "$work"

run() {  # name n_vox iti_min iti_max
  local dir="$work/$1"
  Rscript "$here/simulate.R" "$dir" "$2" 1 "$3" "$4"
  python3 "$here/python_lss.py" "$dir"
  Rscript "$here/run_comparison.R" "$dir" 3
}

run rapid_2k 2000 2 8
run slow_2k 2000 8 14
run rapid_20k 20000 2 8

Rscript -e '
args <- commandArgs(TRUE)
files <- file.path(args[1], c("rapid_2k", "slow_2k", "rapid_20k"), "comparison.csv")
tabs <- lapply(files, function(f) cbind(scenario = basename(dirname(f)), read.csv(f)))
write.csv(do.call(rbind, tabs), file.path(args[1], "all_results.csv"), row.names = FALSE)
' "$work"
