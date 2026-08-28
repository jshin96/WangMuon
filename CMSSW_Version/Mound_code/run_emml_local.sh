#!/bin/bash
set -eo pipefail

# --------------------------- EDIT THIS BLOCK ---------------------------
# The reconstruction configuration is intentionally set in this file.  Edit
# these values, save, then run: ./run_emml_local.sh
INPUT_DIR=Incline_5m_pm30deg_Detector
SKIM_FACTOR=1
MPI_PROCS=8
MODE=reconstruct                 # reconstruct or diagnostics
INITIAL_LAMBDA=1.2e-6            # rad^2/cm
MAX_ACCEPTED_MUONS=0       # after detector acceptance; 0 means no cap

# Reconstruction voxel size along global Y. Tracks in this setup travel mostly
# along +Y, so 2 m Y voxels improve event support while retaining four cells
# across the 8 m room. X and Z remain at 1 m.
VOXEL_SIZE_Y_M=2

# Choose exactly one virtual detector selection.  The configurable array is

# Configurable array: number of opposing pairs and individual panel size.
# A muon may enter through one stack and exit through a different stack.
# This applies virtual-panel acceptance to hits recorded on the continuous
# full-2pi inner and outer detector shells.
USE_CONFIGURABLE_DETECTORS=1 #1:on, 0:off
DETECTOR_PAIRS=2
# One azimuth per opposing detector pair, in degrees from the tilted-frame
# +Y axis. Positive angles rotate toward +X. Each entry automatically places
# the opposing stack at angle + 180 degrees.
DETECTOR_PAIR_ANGLES_DEG=(-10 )
DETECTOR_WIDTH_M=4
DETECTOR_HEIGHT_M=1
# Panel-centre distance along the tilted cylinder axis, measured from its
# lower cap in the ground plane.
DETECTOR_ELEVATION_M=4
# Must equal MOUND_GROUND_INCLINE_DEG used by run_sim_job.sh.  This is the
# terrain inclination and the rigid detector cylinder uses this same angle.
GROUND_INCLINE_DEG=15.0
DETECTOR_INCLINE_DEG=${GROUND_INCLINE_DEG}
# Direction cut applied before detector-panel acceptance. This is the incoming
# horizontal unit-direction Y component, dy/sqrt(dx^2+dy^2), in [-1, 1].
#   0.0  accepts only +Y/non-negative-Y travel
#  -0.05 accepts all +Y travel plus a small component toward -Y
#  -1.0  disables the direction cut
MIN_MUON_HORIZONTAL_Y_COMPONENT=0.0

# diagnostics option
DIAGNOSTIC_LAMBDAS="5e-7 1e-6 1.5e-6"


EXTRA_ARGS=""                    # e.g. "--closure 12345"; reconstruction only
SHOW_OPTIONS=0                   # change to 1 to print the guide and exit
SAMPLES=(AirRoom SoilRoom)


EXTRA_EXECUTABLE_ARGS=()
EXTRA_ARGS_TAG=""
if [[ -n "${EXTRA_ARGS}" ]]; then
    read -r -a EXTRA_EXECUTABLE_ARGS <<< "${EXTRA_ARGS}"
    for extra_arg in "${EXTRA_EXECUTABLE_ARGS[@]}"; do
        extra_arg=${extra_arg//-/}
        EXTRA_ARGS_TAG+="${EXTRA_ARGS_TAG:+_}${extra_arg}"
    done
fi

INCLINE_TAG=${DETECTOR_INCLINE_DEG//-/minus}
INCLINE_TAG=${INCLINE_TAG//./p}
ELEVATION_TAG=${DETECTOR_ELEVATION_M//-/minus}
ELEVATION_TAG=${ELEVATION_TAG//./p}
GROUND_TAG=${GROUND_INCLINE_DEG//-/minus}
GROUND_TAG=${GROUND_TAG//./p}
DIRECTION_TAG=${MIN_MUON_HORIZONTAL_Y_COMPONENT//-/minus}
DIRECTION_TAG=${DIRECTION_TAG//./p}
VOXEL_Y_TAG=${VOXEL_SIZE_Y_M//-/minus}
VOXEL_Y_TAG=${VOXEL_Y_TAG//./p}
PAIR_ANGLES_TAG=""
for pair_angle in "${DETECTOR_PAIR_ANGLES_DEG[@]}"; do
    pair_angle_tag=${pair_angle//-/minus}
    pair_angle_tag=${pair_angle_tag//./p}
    PAIR_ANGLES_TAG+="${PAIR_ANGLES_TAG:+_}${pair_angle_tag}"
done
if [[ "${USE_CONFIGURABLE_DETECTORS}" == "1" ]]; then
    OUTPUT_TAG=Detector${DETECTOR_PAIRS}Pair_AnglesFromY${PAIR_ANGLES_TAG}deg_${DETECTOR_WIDTH_M}x${DETECTOR_HEIGHT_M}m_Elevation${ELEVATION_TAG}m_Incline${INCLINE_TAG}deg_MinYDirection${DIRECTION_TAG}_GroundIncline${GROUND_TAG}deg_VoxelY${VOXEL_Y_TAG}m_Skim${SKIM_FACTOR}_${MODE}${EXTRA_ARGS_TAG:+_${EXTRA_ARGS_TAG}}
else
    OUTPUT_TAG=DetectorFull2Pi_MinYDirection${DIRECTION_TAG}_GroundIncline${GROUND_TAG}deg_VoxelY${VOXEL_Y_TAG}m_Skim${SKIM_FACTOR}_${MODE}${EXTRA_ARGS_TAG:+_${EXTRA_ARGS_TAG}}
fi

print_emml_options() {
    printf '%s\n' \
        'run_emml_local.sh controls ./bin/emml_mpi_NoTV as follows:' \
        '' \
        '  Required executable arguments (set by this script)' \
        '    <input_root_file> [skim_factor]' \
        '    INPUT_DIR            Directory containing the matched ROOT files.' \
        '    SKIM_FACTOR          Process every Nth input event (default: 1).' \
        '' \
        '  Reconstruction-likelihood options' \
        '    --angle-only         Use angular scattering only; enabled by this runner.' \
        '    --initial-lambda L   Initial lambda [rad^2/cm]; set INITIAL_LAMBDA.' \
        '    --max-accepted-muons N  Cap accepted tracks before EMML; set MAX_ACCEPTED_MUONS (0 = no cap).' \
        '    --voxel-size-y M   Global-Y voxel size in metres; set VOXEL_SIZE_Y_M (recommended: 2).' \
        '                       X and Z voxel sizes remain 1 m.' \
        '' \
        '  Virtual detector selections (choose exactly one)' \
        '    none                 Full cylindrical detector coverage.' \
        '    configurable array   Default: any number of opposing hillside pairs with configurable panels.' \
        '                         Set DETECTOR_PAIRS, DETECTOR_WIDTH_M, DETECTOR_HEIGHT_M,' \
        '                         DETECTOR_ELEVATION_M, and DETECTOR_PAIR_ANGLES_DEG.' \
        '                         Supply exactly one angle for every pair.' \
        '                         Angles are measured from tilted-frame +Y; positive rotates toward +X.' \
        '                         Every angle creates two opposing detector stacks.' \
        '                         Elevation is axial distance above the lower cap, in metres.' \
        '                         The cylinder inclination is always GROUND_INCLINE_DEG.' \
        '                         The runner passes --detector-pairs N --detector-width W --detector-height H,' \
        '                         --detector-elevation M --detector-incline-deg D, plus one' \
        '                         --detector-pair-angle-deg A per pair.' \
        '                         An entry stack and exit stack may belong to different pairs.' \
        '    MIN_MUON_HORIZONTAL_Y_COMPONENT' \
        '                         Minimum incoming dy/sqrt(dx^2+dy^2), from -1 to 1.' \
        '                         Use 0 for +Y only, a small negative value to allow slight -Y, or -1 to disable.' \
        '    GROUND_INCLINE_DEG   Terrain angle passed as --ground-incline-deg; must match the simulation.' \
        '                         The voxel Z bounds are derived from this angle.' \
        '' \
        '  Validation / diagnostic options' \
        '    MODE=diagnostics     Trace paths and print residual diagnostics; no reconstruction output.' \
        '                         Set DIAGNOSTIC_LAMBDAS="..." for the uniform reference lambdas.' \
        '    --closure [seed]     Replace measured residuals with a synthetic 10x cavity phantom.' \
        '                         Use only as a closure test, never for the GEANT4 room comparison.' \
        '' \
        '  Runner controls' \
        '    MPI_PROCS            MPI ranks (default: 8).' \
        '    MAX_ACCEPTED_MUONS   Maximum post-acceptance tracks passed to EMML.' \
        '    OUTPUT_TAG           Required unique suffix for ROOT/VTK/log outputs.' \
        '    EXTRA_ARGS           Advanced whitespace-separated reconstruction flags; e.g. "--closure 12345".' \
        '' \
        'Edit the configuration block at the top of this file, then run: ./run_emml_local.sh'
}

if [[ "${SHOW_OPTIONS}" == "1" ]]; then
    print_emml_options
    exit 0
fi

selection_count=0
[[ "${USE_CONFIGURABLE_DETECTORS}" == "1" ]] && ((selection_count += 1))
[[ "${USE_VIRTUAL_FOUR_PATCH}" == "1" ]] && ((selection_count += 1))
[[ "${USE_SPARSE_FOUR_PAIR}" == "1" ]] && ((selection_count += 1))
[[ "${USE_SPARSE_EIGHT_PAIR}" == "1" ]] && ((selection_count += 1))
if (( selection_count > 1 )); then
    echo "Error: choose only one virtual detector selection mode." >&2
    exit 1
fi

if [[ "${USE_CONFIGURABLE_DETECTORS}" == "1" ]]; then
    if ! [[ "${DETECTOR_PAIRS}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: DETECTOR_PAIRS must be a positive integer." >&2
        exit 1
    fi
    if (( ${#DETECTOR_PAIR_ANGLES_DEG[@]} != DETECTOR_PAIRS )); then
        echo "Error: DETECTOR_PAIRS=${DETECTOR_PAIRS} requires exactly ${DETECTOR_PAIRS} entries in DETECTOR_PAIR_ANGLES_DEG, but ${#DETECTOR_PAIR_ANGLES_DEG[@]} were provided." >&2
        exit 1
    fi
fi

# These arguments are shared by reconstruction and diagnostics, ensuring that
# diagnostics measures the same accepted paths used by EMML.
DETECTOR_SELECTION_ARGS=()
if [[ "${USE_CONFIGURABLE_DETECTORS}" == "1" ]]; then
    DETECTOR_SELECTION_ARGS=(--detector-pairs "${DETECTOR_PAIRS}"
                             --detector-width "${DETECTOR_WIDTH_M}"
                             --detector-height "${DETECTOR_HEIGHT_M}"
                             --detector-elevation "${DETECTOR_ELEVATION_M}"
                             --detector-incline-deg "${DETECTOR_INCLINE_DEG}")
    for pair_angle in "${DETECTOR_PAIR_ANGLES_DEG[@]}"; do
        DETECTOR_SELECTION_ARGS+=(--detector-pair-angle-deg "${pair_angle}")
    done
elif [[ "${USE_SPARSE_EIGHT_PAIR}" == "1" ]]; then
    DETECTOR_SELECTION_ARGS=(--sparse-eight-pair)
elif [[ "${USE_SPARSE_FOUR_PAIR}" == "1" ]]; then
    DETECTOR_SELECTION_ARGS=(--sparse-four-pair)
elif [[ "${USE_VIRTUAL_FOUR_PATCH}" == "1" ]]; then
    DETECTOR_SELECTION_ARGS=(--virtual-four-patch)
fi
TERRAIN_ARGS=(--ground-incline-deg "${GROUND_INCLINE_DEG}")
DIRECTION_ARGS=(--min-muon-horizontal-y "${MIN_MUON_HORIZONTAL_Y_COMPONENT}")
VOXEL_ARGS=(--voxel-size-y "${VOXEL_SIZE_Y_M}")

# Keep the comparison output free of CMSSW CUDA/RDMA runtime warnings.  The
# complete unfiltered stream remains available in the corresponding .raw.log.
filter_runtime_warnings() {
    awk '
        /^!!!!!!!!!!!!!!!!/ { skip_nvidia = 1; next }
        skip_nvidia && /^Linked to libnvidia-ml library at wrong path :/ {
            skip_nvidia = 0; next
        }
        skip_nvidia { next }
        /^libibverbs: Warning: couldn.t open config directory/ { next }
        { print }
    ' "$1"
}

echo "Starting EMML job on node: $(hostname) | Job ID: ${JOB_ID}"
echo "EMML skim factor: ${SKIM_FACTOR}"
echo "EMML MPI processes: ${MPI_PROCS}"
echo "EMML mode: ${MODE}"
echo "Maximum post-acceptance muons: ${MAX_ACCEPTED_MUONS}"
echo "Terrain inclination: ${GROUND_INCLINE_DEG} deg"
echo "Minimum incoming horizontal Y direction: ${MIN_MUON_HORIZONTAL_Y_COMPONENT}"
echo "Reconstruction voxel size X/Y/Z: 1 / ${VOXEL_SIZE_Y_M} / 1 m"
if [[ "${MODE}" == "reconstruct" ]]; then
    echo "Angle-only initial lambda: ${INITIAL_LAMBDA} rad^2/cm"
    echo "Virtual four-patch acceptance: ${USE_VIRTUAL_FOUR_PATCH}"
    echo "Sparse four-pair 2 m setup: ${USE_SPARSE_FOUR_PAIR}"
    echo "Sparse eight-pair 2 m setup: ${USE_SPARSE_EIGHT_PAIR}"
    echo "Configurable detector array: ${USE_CONFIGURABLE_DETECTORS}"
    if [[ "${USE_CONFIGURABLE_DETECTORS}" == "1" ]]; then
        echo "Detector pairs / panel width / panel height: ${DETECTOR_PAIRS} / ${DETECTOR_WIDTH_M} m / ${DETECTOR_HEIGHT_M} m"
        echo "Detector pair angles from tilted-frame +Y: ${DETECTOR_PAIR_ANGLES_DEG[*]} deg"
        echo "Detector axial centre / cylinder incline: ${DETECTOR_ELEVATION_M} m / ${DETECTOR_INCLINE_DEG} deg"
        echo "Entry and exit may use any two distinct detector stacks."
    fi
elif [[ "${MODE}" == "diagnostics" ]]; then
    echo "Diagnostic lambda scan: ${DIAGNOSTIC_LAMBDAS} rad^2/cm"
    echo "Diagnostics use the same detector selection as reconstruction."
    if [[ -n "${EXTRA_ARGS}" ]]; then
        echo "Error: EXTRA_ARGS is supported only in reconstruction mode." >&2
        exit 1
    fi
else
    echo "Error: MODE must be reconstruct or diagnostics." >&2
    exit 1
fi

# Source the CVMFS environments
source /cvmfs/cms.cern.ch/cmsset_default.sh

# !!! REMEMBER TO USE YOUR ACTUAL ABSOLUTE PATH HERE !!!
cd ${CMSSW_BASE}/src
eval `scramv1 runtime -sh`

# Return to Condor scratch directory where the input file was transferred
cd ${CMSSW_BASE}/src/WangMuon/CMSSW_Version/bin
BUILD_LOG=../${INPUT_DIR}/EMML_NoTV_build.log
if ! mpicxx run_emml_MPI_NoTV.cc MoundTomographyEMML_MPI_NoTV.cc -o emml_mpi_NoTV \
    `root-config --cflags --libs` -O3 > "${BUILD_LOG}" 2>&1; then
    echo "EMML build failed; see ${BUILD_LOG}" >&2
    exit 1
fi
#mpicxx run_emml_MPI_TV.cc MoundTomographyEMML_MPI_TV.cc -o emml_mpi_TV `root-config --cflags --libs` -O3
cd ..

for sample in "${SAMPLES[@]}"; do
    INPUT_FILE=${INPUT_DIR}/MoundTomographyData_${sample}.root
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "Error: missing input file ${INPUT_FILE}" >&2
        exit 1
    fi

    if [[ "${MODE}" == "reconstruct" ]]; then
        OUTPUT_ROOT=${INPUT_DIR}/EMML_Results_${sample}_${OUTPUT_TAG}.root
        OUTPUT_VTK=${INPUT_DIR}/EMML_Results_${sample}_${OUTPUT_TAG}.vtk
        LOG_FILE=${INPUT_DIR}/EMML_Reconstruction_${sample}_${OUTPUT_TAG}.log
        RAW_LOG=${INPUT_DIR}/EMML_Reconstruction_${sample}_${OUTPUT_TAG}.raw.log
#        if [[ -e "${OUTPUT_ROOT}" || -e "${OUTPUT_VTK}" ]]; then
#            echo "Error: output already exists for ${sample}; choose another OUTPUT_TAG." >&2
#            exit 1
#        fi
        if [[ -e EMML_TomographyResults.root || -e EMML_TomographyResults.vtk ]]; then
            echo "Error: found unmoved EMML_TomographyResults output in the project directory." >&2
            exit 1
        fi
        echo "Reconstructing ${sample} with angle-only EMML"
        RECONSTRUCTION_ARGS=(--angle-only --initial-lambda "${INITIAL_LAMBDA}"
                             --max-accepted-muons "${MAX_ACCEPTED_MUONS}"
                             "${TERRAIN_ARGS[@]}"
                             "${DIRECTION_ARGS[@]}"
                             "${VOXEL_ARGS[@]}"
                             "${DETECTOR_SELECTION_ARGS[@]}")
        if [[ -n "${EXTRA_ARGS}" ]]; then
            RECONSTRUCTION_ARGS+=("${EXTRA_EXECUTABLE_ARGS[@]}")
            echo "Additional emml_mpi_NoTV arguments: ${EXTRA_ARGS}"
        fi
        if ! OMPI_MCA_opal_cuda_support=0 mpiexec -np "${MPI_PROCS}" \
            --mca btl tcp,self \
            --mca btl_openib_allow_ib 0 --mca pml ob1 ./bin/emml_mpi_NoTV \
            "${INPUT_FILE}" "${SKIM_FACTOR}" "${RECONSTRUCTION_ARGS[@]}" \
            > "${RAW_LOG}" 2>&1; then
            filter_runtime_warnings "${RAW_LOG}" | tee "${LOG_FILE}"
            echo "EMML run failed; see ${RAW_LOG}" >&2
            exit 1
        fi
        filter_runtime_warnings "${RAW_LOG}" | tee "${LOG_FILE}"
        rm EMML_TomographyResults.root 
        mv EMML_TomographyResults.vtk "${OUTPUT_VTK}"
    else
        for lambda in ${DIAGNOSTIC_LAMBDAS}; do
            label=${lambda//[^0-9A-Za-z]/_}
            LOG_FILE=${INPUT_DIR}/ResidualDiagnostics_${sample}_${OUTPUT_TAG}_lambda_${label}.log
            RAW_LOG=${INPUT_DIR}/ResidualDiagnostics_${sample}_${OUTPUT_TAG}_lambda_${label}.raw.log
            echo "Running ${sample} at lambda=${lambda} rad^2/cm"
            if ! OMPI_MCA_opal_cuda_support=0 mpiexec -np "${MPI_PROCS}" \
                --mca btl tcp,self \
                --mca btl_openib_allow_ib 0 --mca pml ob1 ./bin/emml_mpi_NoTV \
                "${INPUT_FILE}" "${SKIM_FACTOR}" --diagnostics-lambda "${lambda}" --angle-only \
                --max-accepted-muons "${MAX_ACCEPTED_MUONS}" \
                "${TERRAIN_ARGS[@]}" \
                "${DIRECTION_ARGS[@]}" \
                "${VOXEL_ARGS[@]}" \
                "${DETECTOR_SELECTION_ARGS[@]}" \
                > "${RAW_LOG}" 2>&1; then
                filter_runtime_warnings "${RAW_LOG}" | tee "${LOG_FILE}"
                echo "EMML run failed; see ${RAW_LOG}" >&2
                exit 1
            fi
            filter_runtime_warnings "${RAW_LOG}" | tee "${LOG_FILE}"
        done
    fi
done
rm ${INPUT_DIR}/*.log
echo "EMML Job ${JOB_ID} complete!"
