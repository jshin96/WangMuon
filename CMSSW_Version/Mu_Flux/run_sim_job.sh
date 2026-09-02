#!/bin/bash
# This recipe prepares one simulation job and runs the muon spectrometer.
set -euo pipefail
# These two numbers give every Condor job its own output-file name.
cluster_id=${1:-local}; process_id=${2:-0}; job_id="${cluster_id}_${process_id}"
# Geant4 otherwise starts every independent batch process from the same
# default random state. Condor's (Cluster, Process) pair gives each job a
# distinct, reproducible stream. Keep this below Ranecu's positive seed range.
if [[ "${cluster_id}" =~ ^[0-9]+$ && "${process_id}" =~ ^[0-9]+$ ]]; then
    random_seed=$(( (10#${cluster_id} * 1009 + 10#${process_id}) % 2147483561 + 1 ))
else
    # Convenient for a local test; production Condor IDs always take the path above.
    random_seed=1
fi
export MUON_RANDOM_SEED="${random_seed}"
# Put the CMSSW programs, including muon_flux_spectrometer, on PATH.
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_BASE}/src"; eval "$(scramv1 runtime -sh)"; cd "${_CONDOR_SCRATCH_DIR:-$PWD}"


# Ideal: manually align Beam Z and detector Z using X*tan(90-ZENITH DEG)


# Mu- beam travels in +X; energy is uniformly sampled in this interval (GeV kinetic energy).
# Mu- beam covers MUON_BEAM_SIZE_CM by MUON_BEAM_SIZE_CM square surface in Y-Z plane. 
export MUON_BEAM_MIN_ENERGY_GEV=1  # The smallest energy for a pretend muon.
export MUON_BEAM_MAX_ENERGY_GEV=100
export MUON_BEAM_SOURCE_X_M=-4.5
export MUON_BEAM_SOURCE_Y_M=0.0
export MUON_BEAM_SOURCE_Z_M=0.0
export MUON_BEAM_SIZE_CM=8
export MUON_BEAM_THETA_HALF_SPREAD_DEG=0.0
export MUON_BEAM_PHI_HALF_SPREAD_DEG=0.0
export MUON_BEAM_ZENITH_DEG=90.0



# 10 cm square GEMs with nominal normal +X. In2 is -50 cm and In1 -70 cm by default.
# Change these numbers here when you want to move the boards.
# Per-board alignment controls. Coordinates are global metres. Zenith is the
# angle between the detector normal and global +Z; 90 deg retains the nominal
# board orientation (normal +X). Edit these assigned values for a scan.
export MUON_GEM_SIZE_CM="10"
export MUON_GEM_THICKNESS_MM="3"
export MUON_GEM_IN1_X_M="-4.0"
export MUON_GEM_IN1_Y_M="0"
export MUON_GEM_IN1_Z_M="0"
export MUON_GEM_IN1_ZENITH_DEG="85"
export MUON_GEM_IN2_X_M="-0.50"
export MUON_GEM_IN2_Y_M="0"
export MUON_GEM_IN2_Z_M="0"
export MUON_GEM_IN2_ZENITH_DEG="85"
export MUON_GEM_OUT_X_M="2.00"
export MUON_GEM_OUT_Y_M="0"
export MUON_GEM_OUT_Z_M="0"
export MUON_GEM_OUT_ZENITH_DEG="85"




# Helmholtz central uniform-field approximation: a cylinder centred at
# (0,0,0), with axis +Y and field +Y. APERTURE is the cylinder radius.
export MUON_HELMHOLTZ_FIELD_T=1.5
export MUON_HELMHOLTZ_LENGTH_M=0.50
export MUON_HELMHOLTZ_APERTURE_CM=30
export GEM_INTRINSIC_POSITION_UM=500
export MUON_RECO_MIN_BEND_SIGMA=3.0
echo "Beam +X, ${MUON_BEAM_MIN_ENERGY_GEV}-${MUON_BEAM_MAX_ENERGY_GEV} GeV; GEM X [m]: ${MUON_GEM_IN1_X_M}, ${MUON_GEM_IN2_X_M}, ${MUON_GEM_OUT_X_M}"
echo "GEM alignment (Y m, Z m, zenith deg): In1=(${MUON_GEM_IN1_Y_M}, ${MUON_GEM_IN1_Z_M}, ${MUON_GEM_IN1_ZENITH_DEG}) In2=(${MUON_GEM_IN2_Y_M}, ${MUON_GEM_IN2_Z_M}, ${MUON_GEM_IN2_ZENITH_DEG}) Out=(${MUON_GEM_OUT_Y_M}, ${MUON_GEM_OUT_Z_M}, ${MUON_GEM_OUT_ZENITH_DEG})"
echo "Field: ${MUON_HELMHOLTZ_FIELD_T} T along +Y, length ${MUON_HELMHOLTZ_LENGTH_M} m"
echo "Geant4 random seed: ${MUON_RANDOM_SEED} (Condor Cluster=${cluster_id}, Process=${process_id})"
# Ask Geant4 to follow the muons.  run.mac says how many to make.
muon_flux_spectrometer run.mac
output="MuonBeamSpectrometer_${job_id}.root"; test -s MuonBeamSpectrometer.root; mv MuonBeamSpectrometer.root "${output}"
echo "Job complete: ${output}"
