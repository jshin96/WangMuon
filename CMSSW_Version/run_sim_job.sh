#!/bin/bash
# 1. Catch the Condor Process ID
Cluster_ID=$1
Process_ID=$2
JOB_ID=${Cluster_ID}_${Process_ID}
ROOM_MATERIAL=${3:-air}
ROOM_MATERIAL=$(printf '%s' "${ROOM_MATERIAL}" | tr '[:upper:]' '[:lower:]')
case "${ROOM_MATERIAL}" in
    air)  ROOM_TAG=AirRoom ;;
    soil) ROOM_TAG=SoilRoom ;;
    *)
        echo "Error: room material must be 'air' or 'soil', got '${ROOM_MATERIAL}'." >&2
        exit 2
        ;;
esac
echo "Job started on node: $(hostname) | Job ID: ${JOB_ID}"

# 2. Source the CVMFS environments
source /cvmfs/cms.cern.ch/cmsset_default.sh

# !!! REMEMBER TO USE YOUR ACTUAL PATH HERE !!!
cd ${CMSSW_BASE}/src
eval `scramv1 runtime -sh`

# 3. Return to Condor scratch directory
cd ${_CONDOR_SCRATCH_DIR}

# 4. Run the simulation
# Positive values make the terrain descend as global +Y increases.  Override
# this environment variable with 0 for the horizontal-baseline comparison.
export MOUND_GROUND_INCLINE_DEG=15
export MOUND_INCLUDE_ROOM=1
export MOUND_ROOM_MATERIAL="${ROOM_MATERIAL}"
# Generate only from the uphill source arc capable of reaching the requested
# downhill detector window. Angles use the same convention as EMML: zero is
# tilted-frame +Y and positive rotates toward +X. The axial window can later
# be narrowed to the physical panel height/elevation.
MOUND_DETECTOR_WINDOW_MIN_DEG="-1.1"
MOUND_DETECTOR_WINDOW_MAX_DEG="1.1"
MOUND_DETECTOR_WINDOW_HEIGHT_M="1"
MOUND_DETECTOR_WINDOW_ELEVATION_M="1"
MOUND_MUON_MIN_MOMENTUM_GEV="1"
MOUND_MUON_MAX_MOMENTUM_GEV="100"
#MOUND_DETECTOR_WINDOW_MIN_DEG="-30.0"
#MOUND_DETECTOR_WINDOW_MAX_DEG="30.0"
#MOUND_DETECTOR_WINDOW_HEIGHT_M="5"
#MOUND_DETECTOR_WINDOW_ELEVATION_M="2.5"
#MOUND_MUON_MIN_MOMENTUM_GEV="10"
#MOUND_MUON_MAX_MOMENTUM_GEV="100"
# Bound cheap EcoMug/analytic rejection sampling. An exhausted event is
# counted in the PRIMARY_GENERATION summary and skipped before transport.
MOUND_MAX_GENERATION_ATTEMPTS="10000"

# PrimaryGeneratorAction reads these settings from the child process environment.
export MOUND_DETECTOR_WINDOW_MIN_DEG MOUND_DETECTOR_WINDOW_MAX_DEG
export MOUND_DETECTOR_WINDOW_HEIGHT_M MOUND_DETECTOR_WINDOW_ELEVATION_M
export MOUND_MUON_MIN_MOMENTUM_GEV MOUND_MUON_MAX_MOMENTUM_GEV
export MOUND_MAX_GENERATION_ATTEMPTS
echo "Ground inclination: ${MOUND_GROUND_INCLINE_DEG} deg"
echo "Room fill material: ${MOUND_ROOM_MATERIAL}"
echo "Conditional detector azimuth window: [${MOUND_DETECTOR_WINDOW_MIN_DEG}, ${MOUND_DETECTOR_WINDOW_MAX_DEG}] deg from +Y"
echo "Conditional detector axial elevation / height: ${MOUND_DETECTOR_WINDOW_ELEVATION_M} / ${MOUND_DETECTOR_WINDOW_HEIGHT_M} m"
echo "Muon momentum window: ${MOUND_MUON_MIN_MOMENTUM_GEV}--${MOUND_MUON_MAX_MOMENTUM_GEV} GeV/c"
echo "Maximum EcoMug attempts per beam event: ${MOUND_MAX_GENERATION_ATTEMPTS}"
echo "Firing cosmic rays..."
wang_muon run.mac

# 5. RENAME THE ROOT FILE (The Magic Step)
mv MoundTomographyData.root "MoundTomographyData_${ROOM_TAG}_${JOB_ID}.root"
echo "Job ${JOB_ID} complete!"
