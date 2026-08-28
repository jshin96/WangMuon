#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // -----------------------------------------------------------
    // FILTER 1: Is it a muon?
    // -----------------------------------------------------------
    G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
    if (particle != G4MuonMinus::MuonMinusDefinition() && 
        particle != G4MuonPlus::MuonPlusDefinition()) {
        return; // It's an electron, gamma, etc. Ignore it entirely.
    }
    // Reject any secondary muons created via pair-production or decay
    if (step->GetTrack()->GetParentID() != 0) {
        return;
    }


    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;
    G4String volName = volume->GetName();

    
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();

    // Accumulate actual ionisation loss in the active GEM gas.  This is kept
    // separate from the truth boundary crossing below.
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.0) {
        if (volName == "DetectorIn1") fEventAction->AddEnergyIn1(edep);
        else if (volName == "DetectorIn2") fEventAction->AddEnergyIn2(edep);
        else if (volName == "DetectorOut1") fEventAction->AddEnergyOut1(edep);
        else if (volName == "DetectorOut2") fEventAction->AddEnergyOut2(edep);
        else if (volName == "DetectorOut3") fEventAction->AddEnergyOut3(edep);
    }

    // Truth positions and momenta are recorded only at volume entry.
    if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) return;
    const G4double time = step->GetPreStepPoint()->GetGlobalTime();

    if (volName == "DetectorIn1") fEventAction->SetHitIn1(pos, mom, time);
    else if (volName == "DetectorIn2") fEventAction->SetHitIn2(pos, mom, time);
    else if (volName == "DetectorOut1") fEventAction->SetHitOut1(pos, mom, time);
    else if (volName == "DetectorOut2") fEventAction->SetHitOut2(pos, mom, time);
    else if (volName == "DetectorOut3") fEventAction->SetHitOut3(pos, mom, time);
    else if (volName == "PhysMound") {
        fEventAction->AddHitMound();
    }
    else if (volName == "PhysRoom" || volName == "PhysRockWall") {
        fEventAction->AddHitRoom();
    }
    /*
    if (volName == "DetectorIn1") {
        fEventAction->SetHitIn1(step->GetPreStepPoint()->GetPosition(),
                                step->GetPreStepPoint()->GetMomentum());
    } 
    else if (volName == "DetectorIn2") {
        fEventAction->SetHitIn2(step->GetPreStepPoint()->GetPosition(),
                                step->GetPreStepPoint()->GetMomentum());
    }
    else if (volName == "DetectorOut1") {
        fEventAction->SetHitOut1(step->GetPreStepPoint()->GetPosition(),
                                 step->GetPreStepPoint()->GetMomentum());
    }
    else if (volName == "DetectorOut2") {
        fEventAction->SetHitOut2(step->GetPreStepPoint()->GetPosition(),
                                 step->GetPreStepPoint()->GetMomentum());
    }
    else if (volName == "PhysMound") {
        fEventAction->AddHitMound();
    }
    else if (volName == "PhysRoom" || volName == "PhysRockWall") {
        fEventAction->AddHitRoom();
    }
    */
}
