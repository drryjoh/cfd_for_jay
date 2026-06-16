# cfd_for_jay — Standalone CodeJeNN Additions for detonationFoam

This repository contains the minimum set of files needed to apply CodeJeNN
neural-network transport modifications to a clean checkout of
[detonationFoam](https://github.com/JieSun-pku/detonationFoam).

Nothing here is a fork of detonationFoam. The files in `applications/` are
either drop-in replacements for existing detonationFoam source files or
entirely new files. The two shell scripts automate applying and reversing
those changes.

---

## Prerequisites — Installing OpenFOAM 8

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

---

## Quick Start

```bash
git clone <this-repo> cfd_for_jay
cd cfd_for_jay
./generate_solvers.sh          # clones detonationFoam and patches it
source /path/to/OpenFOAM-8/etc/bashrc
cd detonationFoam/applications/solvers
./Allmake                      # build everything
```

After a successful build, confirm with:

```bash
which generateData
which checkNNPrediction
which detonationFoam_V2.0
```

To undo all modifications and restore the clean clone:

```bash
cd cfd_for_jay
./remove_codejenn.sh
```

---

## Scripts

### `generate_solvers.sh`

```
Usage: ./generate_solvers.sh [DEST_DIR]
       DEST_DIR defaults to ./detonationFoam
```

Performs the following steps in order:

1. **Clone** `git@github.com:JieSun-pku/detonationFoam.git` into `DEST_DIR`
   (skipped if the directory already exists).
2. **Patch** eight modified detonationFoam source files by overwriting them
   with the CodeJeNN-aware versions stored in this repository.
3. **Add** five new NN header files into the solver's `transport/` directory.
4. **Add** three new utility applications (`checkMuNNError`,
   `checkNNPrediction`, `generateData`) into `applications/solvers/`.
5. **Install** the `Allmake` build script at `applications/solvers/Allmake`.

### `remove_codejenn.sh`

```
Usage: ./remove_codejenn.sh [DEST_DIR]
       DEST_DIR defaults to ./detonationFoam
```

Fully reverses `generate_solvers.sh`:

1. Deletes the three utility application directories.
2. Deletes `Allmake`.
3. Removes the five new NN header files from the `transport/` directory.
4. Runs `git checkout --` on all eight modified files, restoring the originals
   from the cloned repository's history.

---

## Repository Layout

```
cfd_for_jay/
├── generate_solvers.sh
├── remove_codejenn.sh
└── applications/
    └── solvers/
        ├── Allmake                              # top-level build script
        │
        ├── checkMuNNError/                      # [NEW] post-processing utility
        │   ├── checkMuNNError.C
        │   └── Make/{files,options}
        │
        ├── checkNNPrediction/                   # [NEW] NN validation utility
        │   ├── checkNNPrediction.C
        │   └── Make/{files,options}
        │
        ├── generateData/                        # [NEW] training data generator
        │   ├── generateData.C
        │   └── Make/{files,options}
        │
        └── detonationFoam_V2.0/
            ├── Make/
            │   └── options                      # [MODIFIED]
            ├── detonationFoam_V2.0.C            # [MODIFIED]
            └── solverTypeNS_mixtureAverage/
                ├── solverTypeNS_mixtureAverage.H  # [MODIFIED]
                ├── rhoYEqn_NS_mixtureAverage.H    # [MODIFIED]
                └── transport/
                    ├── readChemistryProperties.H  # [MODIFIED]
                    ├── createDiffFields.H         # [MODIFIED]
                    ├── createMuFields.H           # [MODIFIED]
                    ├── updateTransProperties.H    # [MODIFIED]
                    ├── codeJeNN_diffusion.H       # [NEW] diffusion NN stub
                    ├── nn_diffusion_interface.H   # [NEW] OpenFOAM↔NN bridge
                    ├── nn_mu_interface.H          # [NEW] OpenFOAM↔NN bridge
                    ├── model_hi_NEW.hpp           # [NEW] viscosity NN (NNHI)
                    └── model_lo_NEW.hpp           # [NEW] viscosity NN (NNLO)
```

---

## What the Modifications Do

All changes are confined to the `NS_mixtureAverage` solver path. The Euler
and `NS_Sutherland` paths are completely untouched. Every NN flag defaults
to `false`, so the solver behaves identically to stock detonationFoam unless
you explicitly enable an NN model in `constant/chemistryProperties`.

### New keys in `constant/chemistryProperties`

| Key | Type | Default | Effect |
|---|---|---|---|
| `useDiffNN` | Switch | `false` | Replace polynomial mixture-averaged diffusion with the NN |
| `useMuNN` | Switch | `false` | Replace Wilke-rule viscosity with the NN |
| `muNNmodel` | word | `NN0` | Select viscosity model: `NNHI` (deeper) or `NNLO` (shallow) |
| `computeMuNNerror` | Switch | `false` | Write per-cell relative error `mu_err` each time-step |

### Modified files — summary of changes

**`Make/options`** — adds `-std=c++14` (required for `constexpr` arrays in
the NN headers).

**`detonationFoam_V2.0.C`** — adds two `#include` directives at file scope,
before `main()`, so the NN interface struct definitions are visible throughout
the compilation unit.

**`solverTypeNS_mixtureAverage.H`** — before the time loop, constructs
`NNDiffInterface` and `NNMuInterface` objects, creates the `activeDiff`
reference (an alias that points to either `Diff` or `DiffNN` depending on
`useDiffNN`), and allocates the `HistoricalErrorCount` diagnostic field when
`useMuNN` is active.

**`rhoYEqn_NS_mixtureAverage.H`** — replaces direct references to `Diff[i]`
with `activeDiff[i]`, allowing the species equations to consume either the
polynomial or NN diffusion coefficients without any branching in the equation
itself.

**`readChemistryProperties.H`** — reads the four new configuration switches
from `chemistryProperties`.

**`createDiffFields.H`** — allocates the `DiffNN` field list (same dimensions
as `Diff`) to hold the NN-predicted diffusion coefficients each time-step.

**`createMuFields.H`** — allocates the `mu_err` scalar field (written with
`AUTO_WRITE`) to record per-cell NN viscosity relative error.

**`updateTransProperties.H`** — the core dispatch logic:
- Skips polynomial binary-diffusion calculation when `useDiffNN` is true.
- Skips Wilke per-species viscosity calculation when `useMuNN` is true
  (unless `computeMuNNerror` is also true, in which case the polynomial is
  computed for error comparison only).
- Calls `nnMuInterface.compute()` or the Wilke loop depending on `useMuNN`.
- Calls `nnDiffInterface.compute()` to populate `DiffNN` when `useDiffNN` is
  true.

### New files — summary

**`codeJeNN_diffusion.H`** — stub diffusion NN for a 9-species H/H2/O/O2/OH/
H2O/HO2/H2O2/Ar system. Replace with CodeJeNN output when a trained model is
available. Currently returns a constant 1 × 10⁻⁵ m²/s placeholder.

**`nn_diffusion_interface.H`** — OpenFOAM struct that wraps the diffusion NN.
Provides `check_species_list()` (called once at startup) and `compute()` (called
each time-step to fill `DiffNN`).

**`nn_mu_interface.H`** — OpenFOAM struct that wraps the two viscosity NNs.
Provides species-order validation and two `compute()` overloads: one that fills
`mixMu` only, and one that also increments `HistoricalErrorCount` for cells
where the raw NN output was non-positive.

**`model_hi_NEW.hpp`** / **`model_lo_NEW.hpp`** — CodeJeNN-generated header-only
viscosity NNs for a H₂/O₂/N₂ system. Input order: `[Y_O2, Y_N2, Y_H2, T]`.
`model_hi_NEW` is a four-layer tanh network (4→12→12→12→1). `model_lo_NEW`
is a two-layer relu/linear network (4→1→1).

### Utility applications

**`generateData`** — generates CSV training data (composition + temperature →
Wilke mixture viscosity) using uniform-Dirichlet sampling optionally blended
with physical mixing-line importance sampling. Reads configuration from
`constant/generateDataProperties`.

**`checkNNPrediction`** — validates both NN viscosity models against the
polynomial Wilke reference over a random sample. Writes a CSV with per-sample
errors and prints mean/max statistics.

**`checkMuNNError`** — post-processing tool; reads `mixMu`, `H2`, `O2`, `N2`,
and `T` from the latest time directory of a completed run and writes
`mu_err_hi` and `mu_err_lo` fields showing cell-by-cell relative NN error
against the solver's own Wilke reference.

---

## Build Order

`Allmake` calls the following targets in sequence:

```
generateData        → wmake
checkNNPrediction   → wmake
checkMuNNError      → wmake
detonationFoam_V2.0 → ./Allwmake   (builds fluxSchemes, DLBFoam, ROUNDSchemes,
                                     dynamicMesh2D, then the solver itself)
```

---

## Appendix — Prompts That Produced This Folder

This folder was created by Claude Code in response to two user prompts.

---

### Prompt 1

> In this directory, there are three folders:
>
> 1. `cfd_for_jay`
>    The target folder where the standalone CodeJeNN additions should be created.
>
> 2. `CodeJeNN_CFD`
>    The modified codebase that contains the CodeJeNN-related changes. This needs to be compared against `detonationFoam`.
>
> 3. `detonationFoam`
>    The reference folder.
>
> `CodeJeNN_CFD` and `detonationFoam` are similar codebases. The goal is to use `cfd_for_jay` to create the bare-minimum standalone CodeJeNN additions needed to modify a clean `detonationFoam_V2.0` checkout.
>
> The issue is that `CodeJeNN_CFD` currently contains both the CodeJeNN modifications and the OpenFOAM-based `detonationFoam` source. We need to keep any source files that differ from `detonationFoam_V2.0` completely separate. Then, we should create a script that merges the standalone CodeJeNN additions from `cfd_for_jay` into a clean clone of `detonationFoam`.
>
> The first tasks are:
>
> 1. Compare:
>
> ```bash
> detonationFoam/applications
> ```
>
> against:
>
> ```bash
> CodeJeNN_CFD/applications
> ```
>
> 2. Create standalone files in `cfd_for_jay` that capture only the minimum required differences introduced by `CodeJeNN_CFD/applications`.
>
> The goal is to separate the CodeJeNN-specific modifications so they can be easily included and compiled into a clean `detonationFoam_V2.0` checkout.
>
> An example workflow, run from inside `cfd_for_jay`, should look like:
>
> ```bash
> git clone git@github.com:JieSun-pku/detonationFoam.git
> ./generate_solvers.sh
> ```
>
> The script `generate_solvers.sh` should:
>
> 1. Clone `detonationFoam`.
> 2. Apply the standalone CodeJeNN solver modifications.
> 3. Add or modify only the minimum files needed to build the CodeJeNN-enabled solvers.
>
> Also implement a second script that reverses this process, removing the CodeJeNN additions and restoring the cloned `detonationFoam` tree to its clean reference state.

---

### Prompt 2

> retry

This was sent after the first response stalled mid-execution. It caused the
assistant to resume and complete the file-writing pass that produced this
repository.

---

### Prompt 3

> Add to the main README.md instructions for how to build all solvers if this
> has not already been added. The README should instruct the user to first
> download and install OpenFOAM 8, then build and install the applications
> provided in this repository.
>
> Copy the tutorials from CodeJeNN_CFD into cfd_for_jay, minimising the files
> as appropriate (removing backups, variants, etc.) so the minimum set of files
> needed to run the tutorials is present.
>
> Create a tutorials/README.md that walks the user through the full CodeJeNN
> workflow: generating training data for the viscosity neural network, verifying
> the network, running a baseline splitter-plate CFD case, and then continuing
> that case with and without the NN model active.
>
> Write a portable submit_splitter_1_2.sh (no hardcoded cluster paths) that
> submits two SLURM jobs continuing from the developed solution at
> t = 3.32×10⁻⁴ s: one baseline run (Wilke viscosity) and one NN-enabled run
> (NNHI model).
>
> Add a summary of this prompt to the main README appendix.

---

### What the assistant did

1. Explored the directory trees of `detonationFoam/applications` and
   `CodeJeNN_CFD/applications` to map their full structure.
2. Ran a recursive diff and read every file that differed or was new in
   `CodeJeNN_CFD`.
3. Identified that the `model_hi_NEW.hpp` and `model_lo_NEW.hpp` NN weights
   lived in `CodeJeNN_CFD/07_cfd_implementation/` and were referenced via an
   extra `-I` flag in `Make/options`; relocated them into the solver's own
   `transport/` directory so the existing include path covers them.
4. Wrote 8 modified files, 5 new NN headers, 3 utility application trees,
   2 shell scripts, and 1 `Allmake` build script into this repository.
