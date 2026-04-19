#!/bin/bash
# submit_new_srs.sh
# Submits array jobs for the three additional spacing ratios added after the
# initial sweep. Does NOT re-run existing SR values.
#
# New SRs: 0.25, 0.325, 0.375
# Total: 11 depths x 3 SRs x 500 trials = 16,500 array tasks (33 jobs)
#
# Output filenames: droplet_output_d<D>_sr0p25_t<T>.csv
#                   droplet_output_d<D>_sr0p325_t<T>.csv
#                   droplet_output_d<D>_sr0p375_t<T>.csv
#
# The C++ must be compiled with the 3-decimal SR filename fix before running.

DEPTHS=(10 15 20 25 30 35 40 45 50 55 60)
NEW_SRS=(0.25 0.325 0.375)

for d in "${DEPTHS[@]}"; do
    for sr in "${NEW_SRS[@]}"; do
        sbatch --export=ALL,DEPTH="$d",SR="$sr" droplet_sim.sh
    done
done
