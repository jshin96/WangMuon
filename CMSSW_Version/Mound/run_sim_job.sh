#!/bin/bash
set -euo pipefail

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

# 4. Simulation configuration
# Geometry: 1 enables the air/soil room; 0 removes the room volume.
export MOUND_INCLUDE_ROOM=1
# Room material: air or soil (set by the third command-line argument).
export MOUND_ROOM_MATERIAL="${ROOM_MATERIAL}"
# Interior room dimensions inside the mound, in metres, along global X/Y/Z.
export MOUND_ROOM_SIZE_X_M="4"
export MOUND_ROOM_SIZE_Y_M="4"
export MOUND_ROOM_SIZE_Z_M="3"
# Surrounding granite/soil wall thickness in cm, applied to every room face.
export MOUND_ROOM_WALL_THICKNESS_CM="50"
# Positive angle means terrain descends as global +Y increases.  Use 0 for a
# horizontal baseline; valid values are strictly between -30 and +30 degrees.
export MOUND_GROUND_INCLINE_DEG=15





# Scattering-plate material (fourth argument): tungsten, lead, copper, iron,
# their chemical symbols, or a Geant4/NIST material name such as G4_Al.
export MOUND_SCATTERING_PLATE_MATERIAL="lead"
# Plate thickness in cm.  The historical variable name is retained even when
# the selected plate material is not tungsten.
export MOUND_TUNGSTEN_THICKNESS_CM="3"


# Flat GEM planes: width is global X; height is global Z; Y is the global
# centre; Z is the bottom-edge elevation above the inclined ground at that Y.
# practically removed IN1. 
export MOUND_GEM_IN1_WIDTH_M="0.000001"
export MOUND_GEM_IN1_HEIGHT_M="0.00000001"
export MOUND_GEM_IN1_Y_M="-27"
export MOUND_GEM_IN1_Z_M="1.0"
export MOUND_GEM_IN2_WIDTH_M="2"
export MOUND_GEM_IN2_HEIGHT_M="2"
export MOUND_GEM_IN2_Y_M="-23"
export MOUND_GEM_IN2_Z_M="1"
export MOUND_GEM_OUT1_WIDTH_M="2"
export MOUND_GEM_OUT1_HEIGHT_M="2"
export MOUND_GEM_OUT1_Y_M="23"
export MOUND_GEM_OUT1_Z_M="0.5"
export MOUND_GEM_OUT2_WIDTH_M="2"
export MOUND_GEM_OUT2_HEIGHT_M="2"
export MOUND_GEM_OUT2_Y_M="25"
export MOUND_GEM_OUT2_Z_M="0.5"
export MOUND_GEM_OUT3_WIDTH_M="2"
export MOUND_GEM_OUT3_HEIGHT_M="2"
export MOUND_GEM_OUT3_Y_M="26.04"
export MOUND_GEM_OUT3_Z_M="0.5"





# EcoMug source surface: a cylinder arc centred on the upstream (-Y) side of
# IN1. Its Z range is elevation above the inclined ground at each sampled
# source position, in metres; the arc half-width is degrees.
export MOUND_SOURCE_CLEARANCE_M="0.1"
export MOUND_SOURCE_Z_MIN_M="0"
export MOUND_SOURCE_Z_MAX_M="4"
export MOUND_SOURCE_UPSTREAM_ARC_HALF_WIDTH_DEG="1.5"
# Direction proposal cone about +Y (phi=90 deg). Exact IN1/IN2 ray checks
# remain the final acceptance, so do not make this cone too narrow.
export MOUND_SOURCE_DIRECTION_PHI_HALF_WIDTH_DEG="3"

# Absolute flux normalization written to RunMetadata.  This is the observed
# local horizontal muon rate in Hz/m^2; 129 is EcoMug's default sea-level
# normalization.  Replace it with a site measurement when available.
export MOUND_HORIZONTAL_MUON_RATE_HZ_M2="129"
# Monte-Carlo points used to integrate the configured source rate.  Raise this
# for a smaller normalization uncertainty, at the cost of startup time.
export MOUND_ECOMUG_RATE_INTEGRATION_POINTS="100000"

# GEM digitisation model.  These settings control the reconstructed GEM
# observables only; Geant4 energy deposition remains the transport truth.
# Reconstruction efficiency before charge threshold (fraction in 0--1).
export GEM_EFFICIENCY="0.96"
# Mean ionisation energy in the Ar/CO2 gas, in eV per primary electron.
export GEM_W_VALUE_EV="26"
# Fano factor controlling primary-ionisation fluctuations.
export GEM_FANO_FACTOR="0.20"
# Mean triple-GEM avalanche gain.
export GEM_MEAN_GAIN="10000"
# Avalanche-gain RMS divided by mean gain (Polya-like width).
export GEM_RELATIVE_GAIN_RMS="0.50"
# Front-end equivalent noise charge, RMS, in electrons.
export GEM_ENC_ELECTRONS="1000"
# Discriminator threshold in fC; hits below it have GEM*_Valid=0.
export GEM_THRESHOLD_FC="0.0001"
# Readout strip/pad pitch in mm.
export GEM_STRIP_PITCH_MM="0.40"
# High-charge, single-coordinate position-resolution limit in micrometres.
export GEM_INTRINSIC_POSITION_UM="1000"
# High-charge timing-resolution limit in ns.
export GEM_INTRINSIC_TIME_NS="4"
# Low-charge discriminator time-walk/jitter scale in ns.
export GEM_THRESHOLD_TIME_JITTER_NS="20"

# Accepted cosmic-muon momentum interval in GeV/c (min must be below max).
export MOUND_MUON_MIN_MOMENTUM_GEV="1"
export MOUND_MUON_MAX_MOMENTUM_GEV="100"
# Maximum EcoMug/analytic sampling attempts per beam event; exhaustion skips
# that event and records it in PRIMARY_GENERATION.
export MOUND_MAX_GENERATION_ATTEMPTS="10000"

# Debug only: set this to any non-empty value to print each accepted primary
# direction.  Leave it unset for production runs.
#export MOUND_PRINT_PRIMARY_DIRECTIONS=1

echo "Ground inclination: ${MOUND_GROUND_INCLINE_DEG} deg"
echo "Room enabled: ${MOUND_INCLUDE_ROOM}; room fill material: ${MOUND_ROOM_MATERIAL}"
echo "Room interior size (X/Y/Z): ${MOUND_ROOM_SIZE_X_M} / ${MOUND_ROOM_SIZE_Y_M} / ${MOUND_ROOM_SIZE_Z_M} m"
echo "Room wall thickness: ${MOUND_ROOM_WALL_THICKNESS_CM} cm"
echo "Scattering-plate material: ${MOUND_SCATTERING_PLATE_MATERIAL}"
echo "Scattering-plate thickness: ${MOUND_TUNGSTEN_THICKNESS_CM} cm"
echo "Flat GEM centres (Y/Z m): In1=${MOUND_GEM_IN1_Y_M}/${MOUND_GEM_IN1_Z_M}, In2=${MOUND_GEM_IN2_Y_M}/${MOUND_GEM_IN2_Z_M}, Out1=${MOUND_GEM_OUT1_Y_M}/${MOUND_GEM_OUT1_Z_M}, Out2=${MOUND_GEM_OUT2_Y_M}/${MOUND_GEM_OUT2_Z_M}, Out3=${MOUND_GEM_OUT3_Y_M}/${MOUND_GEM_OUT3_Z_M}"
echo "EcoMug upstream source: clearance=${MOUND_SOURCE_CLEARANCE_M} m, Z=[${MOUND_SOURCE_Z_MIN_M}, ${MOUND_SOURCE_Z_MAX_M}] m, position/direction half-width=${MOUND_SOURCE_UPSTREAM_ARC_HALF_WIDTH_DEG}/${MOUND_SOURCE_DIRECTION_PHI_HALF_WIDTH_DEG} deg"
echo "Muon flux normalization: ${MOUND_HORIZONTAL_MUON_RATE_HZ_M2} Hz/m^2; EcoMug rate-integration points: ${MOUND_ECOMUG_RATE_INTEGRATION_POINTS}"
echo "GEM response: efficiency=${GEM_EFFICIENCY}, gain=${GEM_MEAN_GAIN}, threshold=${GEM_THRESHOLD_FC} fC, ENC=${GEM_ENC_ELECTRONS} e-"
echo "Muon momentum window: ${MOUND_MUON_MIN_MOMENTUM_GEV}--${MOUND_MUON_MAX_MOMENTUM_GEV} GeV/c"
echo "Maximum EcoMug attempts per beam event: ${MOUND_MAX_GENERATION_ATTEMPTS}"
echo "Firing cosmic rays..."
wang_muon run.mac

# 5. RENAME THE ROOT FILE (The Magic Step)
OUTPUT_FILE="MoundTomographyData_${ROOM_TAG}_${JOB_ID}.root"
if [[ ! -s MoundTomographyData.root ]]; then
    echo "Error: simulation did not create a non-empty MoundTomographyData.root." >&2
    exit 1
fi
mv MoundTomographyData.root "${OUTPUT_FILE}"
echo "Job ${JOB_ID} complete: ${OUTPUT_FILE}"
