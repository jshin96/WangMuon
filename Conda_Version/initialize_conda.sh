#!/usr/bin/env bash
# Create (or update) the Conda environment needed by this WangMuon project.
#
# Usage:
#   ./initialize_conda.sh
#   conda activate WangMuon
#
# This project needs Geant4 (simulation), ROOT (I/O and reconstruction), and
# Open MPI (the EMML reconstruction).  Their compatible runtime dependencies,
# including Geant4 physics data, are resolved by conda-forge.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
environment_file="${script_dir}/environment.yml"
environment_name="WangMuon"

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: Conda was not found on PATH. Install Miniforge or Conda first." >&2
    exit 1
fi

if conda env list | awk 'NR > 2 { print $1 }' | grep -Fxq "${environment_name}"; then
    echo "Updating Conda environment '${environment_name}' from ${environment_file}..."
    conda env update --name "${environment_name}" --file "${environment_file}"
else
    echo "Creating Conda environment '${environment_name}' from ${environment_file}..."
    conda env create --file "${environment_file}"
fi

echo
echo "Environment ready. Start using it with:"
echo "  conda activate ${environment_name}"
echo
echo "Quick checks:"
echo "  geant4-config --version"
echo "  root-config --version"
echo "  mpirun --version"
