# CodeJeNN Tutorials

These tutorials demonstrate the full CodeJeNN workflow: generating training
data for the viscosity neural network, verifying the network, running a
baseline splitter-plate CFD case, and then continuing that case with and
without the NN model active.

---

## Prerequisites — Installing OpenFOAM 8 and Building the Applications

### 1. Install OpenFOAM 8

Download and install OpenFOAM 8 from the official source:

```
https://openfoam.org/download/8-source/
```

Follow the platform-specific build instructions provided there.  After
installation, confirm that sourcing the environment works:

```bash
source /path/to/OpenFOAM-8/etc/bashrc
foamVersion   # should print 8
```

### 2. Clone detonationFoam and apply CodeJeNN modifications

From the root of this repository (`cfd_for_jay/`):

```bash
./generate_solvers.sh
```

This clones `JieSun-pku/detonationFoam` and patches it with the CodeJeNN
modifications.  See the main [README](../README.md) for details.

### 3. Build all applications

```bash
source /path/to/OpenFOAM-8/etc/bashrc
cd detonationFoam/applications/solvers
./Allmake
```

`Allmake` builds the following targets in order:

| Target | Description |
|---|---|
| `generateData` | Training-data generator |
| `checkNNPrediction` | NN validation utility |
| `checkMuNNError` | Post-processing error utility |
| `detonationFoam_V2.0` | Main CFD solver (via upstream `Allwmake`) |

All binaries are installed to `$FOAM_USER_APPBIN`.  After a successful build,
confirm with:

```bash
which generateData
which checkNNPrediction
which detonationFoam_V2.0
```

---

## Tutorial Layout

```
tutorials/
├── submit_splitter_1_2.sh         SLURM submit script (step 5)
├── H2_O2_N2_mu_data/             Step 2 — generate training data
├── H2_O2_N2_check_nn/            Step 3 — verify NN implementation
├── H2_O2_N2_splitter/            Step 4 — run the splitter-plate case
└── H2_O2_N2_splitter_developed/  Pre-computed developed solution (t = 3.32×10⁻⁴ s)
```

All tutorials model a **H₂/air splitter-plate mixing layer** where air is
represented as an N₂/O₂ mixture (Y_N2 = 0.767, Y_O2 = 0.233).  The hydrogen
stream enters at 350 K (M = 0.3) and the air stream at 1500 K (M = 0.3).
Chemistry is disabled; the case is non-reacting.

---

## Step 1 — Source the OpenFOAM environment

All of the following steps require the OpenFOAM environment to be active.
Run this once in each terminal session before proceeding:

```bash
source /path/to/OpenFOAM-8/etc/bashrc
```

---

## Step 2 — Generate NN Training Data

**Directory:** `H2_O2_N2_mu_data/`

This step generates the CSV data used to train the mixture-viscosity neural
networks.  The sampler draws from the physical mixing-line between a pure-H₂
fuel stream and an air oxidizer stream, blended with uniform Dirichlet
sampling to broaden coverage.

```bash
cd H2_O2_N2_mu_data
generateData
```

Output: `mu_training_data.csv` (10 000 rows, columns: `j, H2, O2, N2, T, mu`)

Each row is one (composition, temperature) → mixture viscosity sample computed
from the same log-polynomial Wilke rule used by the solver.  The data is used
to train the viscosity NNs embedded in `model_hi_NEW.hpp` and
`model_lo_NEW.hpp`.

**Sampling parameters** (`constant/generateDataProperties`):

| Parameter | Value | Meaning |
|---|---|---|
| `nSamples` | 10 000 | Total samples |
| `fracPhysical` | 0.4 | Fraction from the H₂/air mixing line |
| `T_fuel` | 350 K | Hydrogen stream temperature |
| `T_oxidizer` | 1500 K | Air stream temperature |
| `T_scatter` | ±200 K | Random scatter around the blended temperature |

---

## Step 3 — Verify the Neural-Network Implementation

**Directory:** `H2_O2_N2_check_nn/`

This step validates both NN viscosity models (NNHI and NNLO) against the
polynomial Wilke reference using 2 000 random samples drawn from the same
distribution as the training data.

```bash
cd H2_O2_N2_check_nn
checkNNPrediction
```

Output: `nn_check_results.csv` and a terminal summary:

```
Results written to 'nn_check_results.csv'
  NNHI | mean: X.XX %  max: X.XX %
  NNLO | mean: X.XX %  max: X.XX %
```

A low mean error (< 1 %) confirms that the NN weights in
`model_hi_NEW.hpp` / `model_lo_NEW.hpp` are correctly embedded and that the
species ordering in the interface (`[Y_O2, Y_N2, Y_H2, T]`) matches what the
solver passes at runtime.

---

## Step 4 — Run the Splitter-Plate Case

**Directory:** `H2_O2_N2_splitter/`

This is the baseline CFD case.  It runs the non-reacting H₂/air splitter-plate
mixing layer for **five flow-through times** (t = 0 → 3.32×10⁻⁴ s, where one
flow-through time is L/U_air = 0.015 m / 226 m/s ≈ 6.64×10⁻⁵ s).

### Mesh and initial field setup

```bash
cd H2_O2_N2_splitter
blockMesh
setFields
```

`blockMesh` generates the base structured mesh from `system/blockMeshDict`.
`setFields` applies the initial species and temperature conditions defined in
`system/setFieldsDict`.

For a refined mesh (matching the developed solution in `H2_O2_N2_splitter_developed`),
use the full preparation script instead:

```bash
./system/prepare_mesh_fields.sh
```

This runs blockMesh followed by three successive passes of topoSet + refineMesh
before calling setFields.

### Decompose and run

```bash
decomposePar
mpirun -n <nprocs> detonationFoam_V2.0 -parallel
```

Replace `<nprocs>` with the value of `numberOfSubdomains` in
`system/decomposeParDict` (default: **128**).

The solver writes one snapshot every flow-through time
(`writeInterval 6.64e-5`).  The run ends at `endTime 3.32e-4`.

### NN switches in `constant/chemistryProperties`

| Key | Default | Description |
|---|---|---|
| `useMuNN` | `off` | Enable NN viscosity (NNHI or NNLO) |
| `muNNmodel` | `NNHI` | Model selection (only active when `useMuNN on`) |
| `useDiffNN` | `off` | Enable NN diffusion (stub — replace with trained model) |
| `computeMuNNerror` | `off` | Write per-cell `mu_err` field each time-step |

The developed solution at t = 3.32×10⁻⁴ s is already provided in
`H2_O2_N2_splitter_developed/0.000332/`.  You can skip step 4 and proceed
directly to step 5 if you only want to run the continuation.

---

## Step 5 — Continue for Five Additional Flow-Through Times (with and without NN)

**Script:** `submit_splitter_1_2.sh`

This script submits two SLURM jobs that both start from the developed solution
at t = 3.32×10⁻⁴ s and run for an additional five flow-through times
(t → 6.64×10⁻⁴ s):

| Job | Directory | `useMuNN` | `muNNmodel` | Description |
|---|---|---|---|---|
| `splitter_noNN` | `H2_O2_N2_splitter_1/` | `off` | — | Wilke polynomial viscosity (baseline) |
| `splitter_NNHI` | `H2_O2_N2_splitter_2/` | `on` | `NNHI` | NN viscosity active |

### Before submitting

Set `OPENFOAM_RC` to the path of your OpenFOAM 8 `etc/bashrc` and ensure the
developed case has been decomposed:

```bash
cd H2_O2_N2_splitter_developed
decomposePar
cd ..
```

Edit the `OPENFOAM_RC` variable at the top of `submit_splitter_1_2.sh` if it
has not already been set in your environment:

```bash
export OPENFOAM_RC=/path/to/OpenFOAM-8/etc/bashrc
```

### Submit

```bash
./submit_splitter_1_2.sh
```

This creates `H2_O2_N2_splitter_1/` and `H2_O2_N2_splitter_2/`, copies the
developed solution into each, edits `constant/chemistryProperties` to set the
correct NN flags, then submits both jobs to the SLURM scheduler.

Each job runs:

```bash
mpirun -np 128 detonationFoam_V2.0 -parallel > log_of 2>&1
```

Solver output is written to `log_of` in each case directory.  SLURM logs go
to `slurm-<jobid>.out` / `.err`.

---

## Summary of Commands

```bash
# Prerequisites (once per session)
source /path/to/OpenFOAM-8/etc/bashrc

# Step 2 — generate training data
cd tutorials/H2_O2_N2_mu_data && generateData

# Step 3 — verify NNs
cd tutorials/H2_O2_N2_check_nn && checkNNPrediction

# Step 4 — run splitter-plate base case
cd tutorials/H2_O2_N2_splitter
blockMesh && setFields && decomposePar
mpirun -n 128 detonationFoam_V2.0 -parallel

# Step 5 — continue with and without NN
cd tutorials
export OPENFOAM_RC=/path/to/OpenFOAM-8/etc/bashrc
./submit_splitter_1_2.sh
```
