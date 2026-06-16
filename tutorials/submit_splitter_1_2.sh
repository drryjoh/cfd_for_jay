#!/bin/bash
# submit_splitter_1_2.sh
#
# Submits two SLURM jobs that continue the splitter-plate simulation for five
# additional flow-through times, one without the NN (baseline) and one with
# the NNHI viscosity model enabled.
#
# Prerequisites:
#   1. OpenFOAM 8 is installed and OPENFOAM_RC points to its etc/bashrc.
#   2. The applications have been built (see the repository README).
#   3. H2_O2_N2_splitter_developed/ has been decomposed:
#         cd H2_O2_N2_splitter_developed && decomposePar
#
# Usage:
#   Set OPENFOAM_RC to the path of your OpenFOAM 8 etc/bashrc, then run:
#       ./submit_splitter_1_2.sh
#
#   Or export it first:
#       export OPENFOAM_RC=/path/to/OpenFOAM-8/etc/bashrc
#       ./submit_splitter_1_2.sh
#
# Output directories created by this script:
#   H2_O2_N2_splitter_1/   — baseline (Wilke viscosity, no NN)
#   H2_O2_N2_splitter_2/   — NN-enabled (NNHI viscosity model)
#
# Each job writes solver output to <CASE_DIR>/log_of and SLURM logs to
# <CASE_DIR>/slurm-<jobid>.out / .err.

# ── Configuration ─────────────────────────────────────────────────────────────

# Directory containing the tutorial cases (default: directory of this script).
BASE="${BASE:-$(cd "$(dirname "$0")" && pwd)}"

# Path to the OpenFOAM 8 etc/bashrc.  Override via the environment variable.
OPENFOAM_RC="${OPENFOAM_RC:-/path/to/OpenFOAM-8/etc/bashrc}"

SRC="${BASE}/H2_O2_N2_splitter_developed"

# ── Helper ────────────────────────────────────────────────────────────────────

submit_case()
{
    local CASE_DIR="$1"
    local FLAG="$2"     # on | off
    local MODEL="$3"    # NNHI | NNLO
    local JOBNAME="$4"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOBNAME}
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --ntasks=128
#SBATCH --output=${CASE_DIR}/slurm-%j.out
#SBATCH --error=${CASE_DIR}/slurm-%j.err

progress() {
    echo ""
    echo "================================================================"
    echo "[\$(date)] \$*"
    echo "================================================================"
}

run_cmd() {
    progress "RUNNING: \$*"
    "\$@"
    local rc=\$?
    progress "FINISHED rc=\${rc}: \$*"
    if [ "\${rc}" -ne 0 ]; then
        progress "FAILED: \$*"
        exit "\${rc}"
    fi
}

progress "Job started"
progress "Job ID: \${SLURM_JOB_ID:-unknown}"
progress "Node list: \${SLURM_JOB_NODELIST:-unknown}"
progress "Case dir: ${CASE_DIR}"
progress "Model: ${MODEL}  useMuNN: ${FLAG}"

progress "Sourcing OpenFOAM"
. ${OPENFOAM_RC}
[ \$? -ne 0 ] && exit 1

progress "Changing into case directory"
run_cmd cd ${CASE_DIR}

progress "Copying developed case into working directory"
run_cmd cp -rv ${SRC}/. .

progress "Editing chemistryProperties  (useMuNN=${FLAG}, muNNmodel=${MODEL})"
sed -i \
  -e 's/useMuNN[[:space:]]\+off;/useMuNN                 ${FLAG};/' \
  -e 's/useMuNN[[:space:]]\+on;/useMuNN                 ${FLAG};/' \
  -e 's/muNNmodel[[:space:]]\+[A-Za-z0-9_]\+;/muNNmodel               ${MODEL};/' \
  constant/chemistryProperties

grep -nE "useMuNN|muNNmodel" constant/chemistryProperties

progress "Starting solver (128 MPI ranks)"
mpirun -np 128 detonationFoam_V2.0 -parallel > log_of 2>&1
solver_rc=\$?

progress "Solver exited rc=\${solver_rc}"
tail -40 log_of || true
[ "\${solver_rc}" -ne 0 ] && exit "\${solver_rc}"
progress "Job completed successfully"
EOF
}

# ── Create output directories and submit ──────────────────────────────────────

mkdir -p "${BASE}/H2_O2_N2_splitter_1"
mkdir -p "${BASE}/H2_O2_N2_splitter_2"

submit_case "${BASE}/H2_O2_N2_splitter_1" "off" "NNHI" "splitter_noNN"
submit_case "${BASE}/H2_O2_N2_splitter_2" "on"  "NNHI" "splitter_NNHI"
