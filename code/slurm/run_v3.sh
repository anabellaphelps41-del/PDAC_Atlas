#!/bin/bash
# run_v3.sh -- submit the v3 rebuild as one dependency chain.
#
#   00m manifest -> 01 harmonize -> 01c doublets -> 02 subset -> 02b batch genes
#                                                              -> 03 TNK | 03 stromal | 03 Bcell
#
# Usage:  bash code/slurm/run_v3.sh            # whole chain
#         bash code/slurm/run_v3.sh 02         # start at step 02 (earlier outputs must exist)
#         AFTER=123 bash code/slurm/run_v3.sh 01c   # start at 01c once job 123 succeeds
#
# Each step only starts if the previous one exited 0 (afterok). Logs go to
# PDAC_Atlas/logs/<jobname>_<jobid>.{out,err}. Steps 04/04b are run by hand.

set -euo pipefail
S=/ix/tbruno/anp795/PDAC_Atlas/code/slurm
FROM=${1:-00m}
order=(00m 01 01c 02 02b 03)
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
    02)  dep=$(submit "$S/02_subset.sbatch") ;;
    02b) dep=$(submit "$S/02b_batchgenes.sbatch") ;;
    03)
      for sub in TNK stromal Bcell; do
        id=$(submit -J "rb03_$sub" "$S/03_pipeline.sbatch" "$sub")
        echo "  03 $sub -> job $id"
      done
      continue ;;
  esac
  echo "  $step -> job $dep"
done

if (( ! started )); then
  echo "unknown start step '$FROM' (expected one of: ${order[*]})" >&2
  exit 1
fi
