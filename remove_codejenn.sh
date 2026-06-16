#!/usr/bin/env bash
# remove_codejenn.sh
#
# Reverses the CodeJeNN modifications applied by generate_solvers.sh.
# Restores the detonationFoam clone to its original clean checkout state.
#
# Run from:  cfd_for_jay/
# Usage:     ./remove_codejenn.sh [DEST_DIR]
#            DEST_DIR defaults to ./detonationFoam

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-${SCRIPT_DIR}/detonationFoam}"

if [ ! -d "${DEST}/.git" ]; then
    echo "ERROR: ${DEST} is not a git repository." >&2
    exit 1
fi

DEST_SOL="${DEST}/applications/solvers"

# ── 1. Remove new utility applications ───────────────────────────────────────
echo "Removing utility applications ..."
for util in checkMuNNError checkNNPrediction generateData; do
    rm -rf "${DEST_SOL}/${util}"
done

rm -f "${DEST_SOL}/Allmake"

# ── 2. Remove new NN files added to the transport directory ──────────────────
echo "Removing NN interface headers ..."
TRANSPORT="${DEST_SOL}/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport"

for f in \
    codeJeNN_diffusion.H \
    nn_diffusion_interface.H \
    nn_mu_interface.H \
    model_hi_NEW.hpp \
    model_lo_NEW.hpp
do
    rm -f "${TRANSPORT}/${f}"
done

# ── 3. Restore modified detonationFoam_V2.0 files from git ───────────────────
echo "Restoring original detonationFoam files from git ..."
(
    cd "${DEST}"
    git checkout -- \
        applications/solvers/detonationFoam_V2.0/Make/options \
        applications/solvers/detonationFoam_V2.0/detonationFoam_V2.0.C \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/solverTypeNS_mixtureAverage.H \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/rhoYEqn_NS_mixtureAverage.H \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createDiffFields.H \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/createMuFields.H \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/readChemistryProperties.H \
        applications/solvers/detonationFoam_V2.0/solverTypeNS_mixtureAverage/transport/updateTransProperties.H
)

echo ""
echo "Clean restore complete: ${DEST}"
echo "The detonationFoam clone is back to its original checkout state."
