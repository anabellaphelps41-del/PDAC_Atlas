#!/bin/bash
# run_v3.sh -- submit the atlas-level steps as one dependency chain.
#
#   00m manifest -> 01 harmonize -> 01c doublets      (outputs in rebuild_v3)
#
# Steps 02 onward are v4 now (DecontX, routing, five subsets): use run_v4.sh.
# The v3 pass of 02-04b is frozen in rebuild_v3 and is not rerun from here.
#
# Usage:  bash code/slurm/run_v3.sh            # whole chain
#         bash code/slurm/run_v3.sh 01c        # start at step 01c (earlier outputs must exist)
#         AFTER=123 bash code/slurm/run_v3.sh 01c   # start at 01c once job 123 succeeds
#
# Each step only starts if the previous one exited 0 (afterok). Logs go to
# PDAC_Atlas/logs/<jobname>_<jobid>.{out,err}.

set -euo pipefail
S=/ix/tbruno/anp795/PDAC_Atlas/code/slurm
FROM=${1:-00m}
order=(00m 01 01c)
started=0
dep="${AFTER:-}"   # AFTER=<jobid> chains the first step onto a job already queued

submit() {  # submit <sbatch args...>; echoes the job id
  local out
  out=$(sbatch --parsable ${dep:+--dependency=afterok:$dep} "$@")
  echo "${out%%;*}"
}

for step in "${order[@]}"; do
  [[ $step == "$FROM" ]] && started=1
  (( started )) || continue
  case $step in
    00m) dep=$(submit "$S/00m_manifest.sbatch") ;;
    01)  dep=$(submit "$S/01_harmonize.sbatch") ;;
    01c) dep=$(submit "$S/01c_doublets.sbatch") ;;
  esac
  echo "  $step -> job $dep"
done

if (( ! started )); then
  echo "unknown start step '$FROM' (expected one of: ${order[*]}; 02 onward: run_v4.sh)" >&2
  exit 1
fi
