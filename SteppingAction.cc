#include "SteppingAction.hh"
#include "G4Step.hh"

SteppingAction::SteppingAction(RunAction* runAction) : fRunAction(runAction) {}
void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    
    // We only care about the primary muon, ignore secondary electrons/photons
    if (track->GetTrackID() != 1) return;

    auto prePoint = step->GetPreStepPoint();
    auto postPoint = step->GetPostStepPoint();
    
    G4String preVol = prePoint->GetPhysicalVolume()->GetName();
    // PostStep volume can be null if particle leaves the world
    G4String postVol = postPoint->GetPhysicalVolume() ? postPoint->GetPhysicalVolume()->GetName() : "OutOfWorld";

    // 1. Check for ENTRY into the mound
    if (preVol == "PhysWorld" && postVol == "PhysMound") {
        fRunAction->fEntryPosX = postPoint->GetPosition().x();
        fRunAction->fEntryPosY = postPoint->GetPosition().y();
        fRunAction->fEntryPosZ = postPoint->GetPosition().z();
        fRunAction->fEntryPx = postPoint->GetMomentum().x();
        fRunAction->fEntryPy = postPoint->GetMomentum().y();
        fRunAction->fEntryPz = postPoint->GetMomentum().z();
        fRunAction->fEntryE = postPoint->GetTotalEnergy();
    }

    // 2. Check for EXIT from the mound back into the world
    if ((preVol == "PhysMound" || preVol == "PhysRoom") && postVol == "PhysWorld") {
        fRunAction->fExitPosX = postPoint->GetPosition().x();
        fRunAction->fExitPosY = postPoint->GetPosition().y();
        fRunAction->fExitPosZ = postPoint->GetPosition().z();
        fRunAction->fExitPx = postPoint->GetMomentum().x();
        fRunAction->fExitPy = postPoint->GetMomentum().y();
        fRunAction->fExitPz = postPoint->GetMomentum().z();
        fRunAction->fExitE = postPoint->GetTotalEnergy();
    }

    // 3. Track Distances
    if (preVol == "PhysMound") {
        fRunAction->fDistDirt += step->GetStepLength();
    } else if (preVol == "PhysRoom") {
        fRunAction->fPassedRoom = 1;
        fRunAction->fDistRoom += step->GetStepLength();
    }

    // 4. Record position & 4-momentum whenever energy changes
    double preE = prePoint->GetKineticEnergy();
    double postE = postPoint->GetKineticEnergy();
    
    if (preE != postE) {
        fRunAction->fStepPosX.push_back(postPoint->GetPosition().x());
        fRunAction->fStepPosY.push_back(postPoint->GetPosition().y());
        fRunAction->fStepPosZ.push_back(postPoint->GetPosition().z());
        fRunAction->fStepPx.push_back(postPoint->GetMomentum().x());
        fRunAction->fStepPy.push_back(postPoint->GetMomentum().y());
        fRunAction->fStepPz.push_back(postPoint->GetMomentum().z());
        fRunAction->fStepE.push_back(postPoint->GetTotalEnergy());
    }
}
