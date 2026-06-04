#!/bin/bash
# 1. Catch the Condor Process ID
Cluster_ID=$1
Process_ID=$2
JOB_ID=${Cluster_ID}_${Process_ID}
echo "Job started on node: $(hostname) | Job ID: ${JOB_ID}"

# 2. Source the CVMFS environments
source /cvmfs/cms.cern.ch/cmsset_default.sh

# !!! REMEMBER TO USE YOUR ACTUAL PATH HERE !!!
cd ${CMSSW_BASE}/src
eval `scramv1 runtime -sh`

# 3. Return to Condor scratch directory
cd ${_CONDOR_SCRATCH_DIR}

# 4. Run the simulation
echo "Firing cosmic rays..."
wang_muon run.mac

# 5. RENAME THE ROOT FILE (The Magic Step)
mv wang_muon_data.root wang_muon_data_${JOB_ID}.root

echo "Job ${JOB_ID} complete!"
