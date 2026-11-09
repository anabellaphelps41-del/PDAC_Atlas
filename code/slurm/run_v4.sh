#!/bin/bash
# run_v4.sh -- submit the v4 pass as one dependency chain.
#
#   02a routing -> 01d DecontX + QC floor -> 02 subsets (5) -> 02b batch genes
#     -> for each subset: 03 Harmony/clusters -> 04 markers -> 04b SingleR gate
#
# Reads 00m/01/01c outputs from rebuild_v3 and writes everything to rebuild_v4
# (see code/rebuild/io.R). 04 and 04b both re-save the subset's harmony object,
# so within a subset they run one after the other, never at the same time.
#
# Usage:  bash code/slurm/run_v4.sh             # whole chain
#         bash code/slurm/run_v4.sh 02          # start at step 02 (earlier outputs must exist)
#         RES=0.8 bash code/slurm/run_v4.sh 04  # 04 -> 04b only, at another resolution
#         AFTER=123 bash code/slurm/run_v4.sh 02b   # start once job 123 succeeds
#         SUBSETS="stromal Bcell" bash code/slurm/run_v4.sh 03
#
# Each step only starts if the previous one exited 0 (afterok). Logs go to
# PDAC_Atlas/logs/<jobname>_<jobid>.{out,err}.

set -euo pipefail
S=/ix/tbruno/anp795/PDAC_Atlas/code/slurm
FROM=${1:-02a}
RES=${RES:-0.6}
SUBSETS=${SUBSETS:-"TNK stromal Bcell myeloid endothelial"}
order=(02a 01d 02 02b 03 04)
started=0
dep="${AFTER:-}"

submit() {  # submit <sbatch args...>; echoes the job id
  local out
  out=$(sbatch --parsable ${dep:+--dependency=afterok:$dep} "$@")
  echo "${out%%;*}"
}

for step in "${order[@]}"; do
  [[ $step == "$FROM" ]] && started=1
  (( started )) || continue
  case $step in
    02a) dep=$(submit "$S/02a_routing.sbatch") ;;
    01d) dep=$(submit "$S/01d_decontx.sbatch") ;;
    02)  dep=$(submit "$S/02_subset.sbatch") ;;
    02b) dep=$(submit -J rb4_02b_batchgenes "$S/02b_batchgenes.sbatch") ;;
    03|04)
      base=$dep
      for sub in $SUBSETS; do
        dep=$base
        if [[ $step == 03 ]]; then
          dep=$(submit -J "rb4_03_$sub" "$S/03_pipeline.sbatch" "$sub")
          echo "  03 $sub -> job $dep"
        fi
        dep=$(submit -J "rb4_04_${sub}_$RES" "$S/04_markers.sbatch" "$sub" "$RES")
        echo "  04 $sub res $RES -> job $dep"
        dep=$(submit -J "rb4_04b_${sub}_$RES" "$S/04b_singler.sbatch" "$sub" "$RES")
        echo "  04b $sub res $RES -> job $dep"
      done
      break ;;
  esac
  echo "  $step -> job $dep"
done

if (( ! started )); then
  echo "unknown start step '$FROM' (expected one of: ${order[*]})" >&2
  exit 1
fi
