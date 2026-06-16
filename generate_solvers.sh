#!/usr/bin/env bash
# generate_solvers.sh
#
# Clones a clean detonationFoam repository and applies the minimum
# CodeJeNN-specific modifications needed to build the NN-enabled solvers.
#
# Run from:  cfd_for_jay/
# Usage:     ./generate_solvers.sh [DEST_DIR]
#            DEST_DIR defaults to ./detonationFoam
#
# After completion, source your OpenFOAM environment and then run:
#   cd DEST_DIR/applications/solvers && ./Allmake
# or to build just the main solver via the upstream Allwmake:
#   cd DEST_DIR/applications/solvers/detonationFoam_V2.0 && ./Allwmake

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-${SCRIPT_DIR}/detonationFoam}"
REPO="git@github.com:JieSun-pku/detonationFoam.git"

SRC="${SCRIPT_DIR}/applications/solvers"
DEST_SOL="${DEST}/applications/solvers"

# ── 1. Clone ──────────────────────────────────────────────────────────────────
if [ -d "${DEST}/.git" ]; then
    echo "INFO: ${DEST} already exists — skipping clone."
else
    echo "Cloning ${REPO} into ${DEST} ..."
    git clone "${REPO}" "${DEST}"
fi

# ── 2. Patch detonationFoam_V2.0 — modified files ────────────────────────────
echo "Patching detonationFoam_V2.0 ..."

# Build system: adds -std=c++14 flag (required for constexpr arrays in NN headers)
cp "${SRC}/detonationFoam_V2.0/Make/options" \
   "${DEST_SOL}/detonationFoam_V2.0/Make/options"

# Main entry point: adds file-scope NN interface includes
cp "${SRC}/detonationFoam_V2.0/detonationFoam_V2.0.C" \
   "${DEST_SOL}/detonationFoam_V2.0/detonationFoam_V2.0.C"

# Solver loop: initialises NN interfaces and activeDiff alias before the time loop
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/solverTypeNS_mixtureAverage.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/solverTypeNS_mixtureAverage.H"

# Species transport: uses activeDiff alias instead of Diff directly
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/rhoYEqn_NS_mixtureAverage.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/rhoYEqn_NS_mixtureAverage.H"

# Transport field creation: adds DiffNN storage for NN-predicted diffusion
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createDiffFields.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createDiffFields.H"

# Mu field creation: adds mu_err diagnostic field
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createMuFields.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createMuFields.H"

# Configuration: reads useDiffNN, useMuNN, muNNmodel, computeMuNNerror from chemistryProperties
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/readChemistryProperties.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/readChemistryProperties.H"

# Transport update: dispatches to NN or polynomial path based on useDiffNN/useMuNN
cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/updateTransProperties.H" \
   "${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/updateTransProperties.H"

# ── 3. Add new CodeJeNN NN headers ────────────────────────────────────────────
echo "Adding NN interface headers ..."
TRANSPORT_DEST="${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport"

for f in \
    codeJeNN_diffusion.H \
    nn_diffusion_interface.H \
    nn_mu_interface.H \
    model_hi_NEW.hpp \
    model_lo_NEW.hpp
do
    cp "${SRC}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/${f}" \
       "${TRANSPORT_DEST}/${f}"
done

# ── 4. Add utility applications ───────────────────────────────────────────────
echo "Adding utility applications ..."
for util in checkMuNNError checkNNPrediction generateData; do
    cp -r "${SRC}/${util}" "${DEST_SOL}/"
done

# ── 5. Install the Allmake build script ───────────────────────────────────────
cp "${SRC}/Allmake" "${DEST_SOL}/Allmake"
chmod +x "${DEST_SOL}/Allmake"

echo ""
echo "CodeJeNN modifications applied to: ${DEST}"
echo ""
echo "To build (OpenFOAM environment must be sourced first):"
echo "   cd ${DEST_SOL} && ./Allmake"
echo ""
echo "To build only the main solver via the upstream Allwmake:"
echo "   cd ${DEST_SOL}/detonationFoam_V2.0 && ./Allwmake"
echo ""
echo "To undo all modifications and restore the clean clone:"
echo "   ./remove_codejenn.sh ${DEST}"
