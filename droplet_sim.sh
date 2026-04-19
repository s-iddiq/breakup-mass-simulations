#!/bin/bash
#SBATCH --job-name=dpm_breakup
#SBATCH --array=0-499
#SBATCH --time=12:00:00
#SBATCH --mem=200M
#SBATCH --partition=pi_ohern
#SBATCH --output=logs/cpp_dbk_d${DEPTH}_sr${SR_TAG}_%A_%a.log
#
# Usage:
#   sbatch --export=ALL,DEPTH=60,SR=0.40 --array=0-499 droplet_sim.sh
#
# DEPTH : obstacle-array height in units of D_obsts  (integer)
# SR    : spacing ratio (min obstacle gap / D_drops)  (float)
# array : trial indices; --array=0-499 gives 500 independent trials

set -uo pipefail

SR_TAG=$(printf "%.2f" "$SR" | tr '.' 'p')
# Strip 3rd decimal zero for new SR values (0.325 stays 0p325, 0.400 -> 0p40)
SR_TAG_3=$(printf "%.3f" "$SR" | tr '.' 'p')
LAST=${SR_TAG_3: -1}
if [ "$LAST" = "0" ]; then
    SR_TAG="${SR_TAG_3%?}"
else
    SR_TAG="$SR_TAG_3"
fi

outdir="results/d${DEPTH}_sr${SR_TAG}"
mkdir -p "$outdir" logs

cd "$outdir"
"$SLURM_SUBMIT_DIR/droplet_sim" "$SLURM_ARRAY_TASK_ID" "$DEPTH" "$SR"
