#!/usr/bin/env bash
# Submit one self-contained Condor simulation campaign.
# Think of this as putting a box of identical simulation jobs into a queue.
#
# Usage:
#   ./submit_sim_run.sh MyRunName
#
# The name becomes a directory beside this script.  Condor returns the ROOT
# output files and logs to that directory, while the copies retained there
# record exactly which job script and macro were submitted.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 RUN_NAME" >&2
    exit 2
fi

RUN_NAME="$1"
if [[ ! "${RUN_NAME}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "Error: RUN_NAME must use only letters, digits, '.', '_' or '-', and must not start with '.'." >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${SCRIPT_DIR}/${RUN_NAME}"

if [[ -e "${RUN_DIR}" ]]; then
    echo "Error: output directory already exists: ${RUN_DIR}" >&2
    echo "Choose a new RUN_NAME so an existing campaign cannot be overwritten." >&2
    exit 1
fi

if ! command -v condor_submit >/dev/null 2>&1; then
    echo "Error: condor_submit is not available in PATH." >&2
    exit 1
fi

mkdir -p "${RUN_DIR}/logs"
cp "${SCRIPT_DIR}/run_sim_job.sh" "${RUN_DIR}/run_sim_job.sh"
cp "${SCRIPT_DIR}/run.mac" "${RUN_DIR}/run.mac"
cp "${SCRIPT_DIR}/submit_sim.sub" "${RUN_DIR}/submit_sim.sub"
cp "${SCRIPT_DIR}/hadder.sh" "${RUN_DIR}/hadder.sh"

# initialdir is both the location for Condor logs and the destination to which
# ON_EXIT file transfer returns the ROOT output.  Override the executable and
# macro input so the exact snapshots above, rather than later edited copies,
# are what the job runs.
(
    cd "${SCRIPT_DIR}"
    condor_submit \
        -append "initialdir = ${RUN_DIR}" \
        -append "executable = ${RUN_DIR}/run_sim_job.sh" \
        -append "transfer_input_files = ${RUN_DIR}/run.mac" \
        -append "+JobBatchName = \"${RUN_NAME}\"" \
        "${SCRIPT_DIR}/submit_sim.sub"
)

echo "Submitted campaign: ${RUN_NAME}"
echo "Results and logs will be returned to: ${RUN_DIR}"
