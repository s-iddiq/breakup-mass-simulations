#!/bin/bash
# submit_all.sh
# Submits one SLURM array job per (depth, SR) combination.
# Total: 11 depths x 5 SRs x 500 trials = 27,500 array tasks
#
# Usage: bash submit_all.sh
#
# To split across McCleary and Grace (they share GPFS):
#   McCleary: --array=0-299
#   Grace:    --array=300-799
# They write to the same results/ directories without collision.

DEPTHS=(10 15 20 25 30 35 40 45 50 55 60)
SRS=(0.20 0.30 0.35 0.40 0.45)

for d in "${DEPTHS[@]}"; do
    for sr in "${SRS[@]}"; do
        sbatch --export=ALL,DEPTH="$d",SR="$sr" droplet_sim.sh
    done
done
