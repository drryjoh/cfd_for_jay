# paper_data — Minimal OpenFOAM Cases for VTK Export and Post-Processing

This directory contains the minimum OpenFOAM case files needed to convert the
final simulation snapshot (t = 6.64×10⁻⁴ s) from both continuation runs into
VTK format and then compare them with `post_error_cfd.py`.

```
paper_data/
├── nn_off/   — baseline case (Wilke polynomial viscosity,  splitter_1)
└── nn_on/    — NN case       (NNHI neural-network viscosity, splitter_2)
```

Each case contains only:

```
<case>/
├── 0.000664/           field snapshots (all fields at the final time step)
├── constant/
│   └── polyMesh/       mesh topology (points, faces, owner, neighbour, boundary)
└── system/
    ├── controlDict     minimal post-processing controlDict (no solver libs)
    ├── fvSchemes
    └── fvSolution
```

---

## Workflow

### 1. Source the OpenFOAM environment

```bash
source /path/to/OpenFOAM-8/etc/bashrc
```

### 2. Convert each case to VTK

Run `foamToVTK` from inside each case directory.  The `-time` flag restricts
conversion to the single stored time step:

```bash
cd paper_data/nn_off
foamToVTK -time 0.000664
cd ../nn_on
foamToVTK -time 0.000664
```

Output is written to `VTK/` inside each case directory.

### 3. Rename the VTK files

`foamToVTK` names its output after the case directory.  Rename the generated
internal-mesh VTK file to match what `post_error_cfd.py` expects:

```bash
# From the paper_data directory:
cp nn_off/VTK/nn_off/nn_off_0.vtk  ../H2_O2_N2_splitter_1.vtk
cp nn_on/VTK/nn_on/nn_on_0.vtk     ../H2_O2_N2_splitter_2.vtk
```

> **Note:** The exact VTK sub-path produced by `foamToVTK` may vary with the
> OpenFOAM version.  Check the `VTK/` directory after conversion and copy
> whichever `.vtk` file represents the internal mesh.

### 4. Run the post-processing script

From the repository root:

```bash
cd /path/to/compare_detonation_foam_cfd_for_jay
python post_error_cfd.py
```

The script reads `H2_O2_N2_splitter_1.vtk` (nn_off baseline) and
`H2_O2_N2_splitter_2.vtk` (nn_on), computes the pointwise relative error in
the x-velocity field, and saves `Ux_error_pct_NN1.pdf` / `.png`.

---

## Case mapping

| Directory  | Source run   | `useMuNN` | VTK filename              |
|------------|--------------|-----------|---------------------------|
| `nn_off/`  | `splitter_1` | `off`     | `H2_O2_N2_splitter_1.vtk` |
| `nn_on/`   | `splitter_2` | `on` (NNHI) | `H2_O2_N2_splitter_2.vtk` |
