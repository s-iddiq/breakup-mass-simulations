# aprilbreakup

DPM (Deformable Particle Method) simulation of droplet breakup through random
obstacle arrays. Parameter sweep over obstacle-array depth and spacing ratio
(SR) on Yale's McCleary HPC cluster.

---

## Project structure

```
aprilbreakup/
├── droplet_sim.cpp          # C++ DPM simulator (compile on cluster)
├── droplet_sim.sh           # SLURM array job script
├── submit_all.sh            # Submit full (depth × SR) sweep
├── submit_new_srs.sh        # Submit additional SR values (0.25, 0.325, 0.375)
├── Analysis.py              # Parse results, filter timeouts, plot histograms
├── fit_beta_delta_zi.py     # Zero-inflated Beta model fitting
└── results/                 # Created at runtime — output CSVs per (depth, SR)
    └── d<D>_sr<SR>/
        └── droplet_output_d<D>_sr<SR>_t<T>.csv
```

---

## Setup

```bash
module load miniconda
conda activate myenv          # or: source ~/myenv/bin/activate
pip install scipy matplotlib numpy --break-system-packages
```

Compile the simulator:

```bash
g++ -O3 -o droplet_sim droplet_sim.cpp -lm
```

---

## Running a sweep

Full sweep (11 depths × 5 SRs × 500 trials = 27,500 jobs):

```bash
bash submit_all.sh
```

Additional SR values:

```bash
bash submit_new_srs.sh
```

Monitor jobs:

```bash
squeue -u sm3546 --format="%.10i %.8T %.5D" | sort
sacct -u sm3546 -S today --format=JobID,State,Elapsed,ExitCode | tail -20
```

**Do not use `scancel -u sm3546`** — this cancels your OOD desktop session.
Cancel specific job IDs instead:

```bash
scancel <jobid1> <jobid2>
```

---

## Parameters

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `depth_cells` | D | Obstacle array height in units of D_obsts |
| `spacing_ratio` | SR | Min gap between obstacles / D_drops |

Current sweep grid:
- Depths: 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60
- SRs: 0.20, 0.25, 0.30, 0.325, 0.35, 0.375, 0.40, 0.45

---

## Analysis

After jobs complete:

```bash
# Parse all results, produce grid histogram + mean vs SR plot
python3 Analysis.py

# Fit zero-inflated Beta model to steady-state data (depths >= 50)
python3 fit_beta_delta_zi.py
```

Outputs:
- `daughter_area_grid.png` — histogram grid (rows=depth, cols=SR)
- `daughter_area_means.png` — mean ± std vs SR per depth
- `daughter_areas_by_depth_sr.csv` — long-format data for downstream use
- `zi_beta_fits.png` — ZI-Beta overlaid on histograms per SR
- `zi_beta_params.png` — w(SR), α(SR), β(SR) parameter trends
- `zi_beta_params.csv` — fitted parameter table

---

## Statistical model

Daughter droplet area distributions are fit with a **zero-inflated Beta**:

```
P(X = 1)  = w                  # intact pass-throughs
X | X < 1 ~ Beta(alpha, beta)  # fragment distribution
```

The delta mass at x=1 captures droplets that traversed the array without
breaking. Area is conserved exactly in the C++ when no breakup occurs, so
intact droplets land at x=1 within floating-point precision. The split
threshold is numerical tolerance (eps=1e-3), not a physical choice.

Parameters are estimated via **method of moments** on the bulk distribution.
MLE is systematically biased by the hard left cutoff at a/A₀ ≈ 0.05, imposed
by the `MIN_VERTICES` floor in the C++ breakup code.

---

## Key findings

- `w` (pass-through fraction) increases monotonically with SR
- Both `alpha` and `beta` decrease with SR — wider gaps produce broader,
  more stochastic fragment distributions
- A gap in a/A₀ ∈ (0.85, 0.99) is a signature of **cascade breakup** dynamics:
  droplets either pass through intact or undergo successive pinching events
  that drive fragments well below 0.85. Single near-intact fragments in this
  range are kinematically rare.
- A transition in concentration parameter α+β around SR ≈ 0.35–0.40 may
  correspond to a crossover between dense-packing and dilute-packing regimes.

---

## Cluster paths

```
Primary:   /gpfs/gibbs/pi/ohern/sm3546/aprilbreakup/
Scratch:   /vast/palmer/scratch/ohern/sm3546/    (purges after 60 days)
```

McCleary and Grace share GPFS — jobs from either cluster write to the same
`results/` directory. Use non-overlapping `--array` ranges to avoid index
collisions (e.g., 0–299 McCleary, 300–799 Grace).
