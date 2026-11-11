#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>

namespace {
G4bool ReadNeutronMode() {
    const char* setting = std::getenv("MOUND_PRIMARY_MODE");
    if (setting == nullptr) return false;
    std::string value(setting);
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value == "neutron";
}
}

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction),
  fNeutronMode(ReadNeutronMode())
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // -----------------------------------------------------------
    // Select the requested particle. Neutron mode also accepts neutron
    // descendants from scattering/cascades, because they contribute to the
    // physically transmitted neutron signal at a downstream detector.
    // -----------------------------------------------------------
    G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
    if (fNeutronMode) {
        if (particle != G4Neutron::NeutronDefinition()) return;
    } else if (particle != G4MuonMinus::MuonMinusDefinition()
               && particle != G4MuonPlus::MuonPlusDefinition()) {
        return;
    }
    // Preserve the existing cosmic-muon definition: only primary muons count.
    if (!fNeutronMode && step->GetTrack()->GetParentID() != 0) {
        return;
    }


    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;
    G4String volName = volume->GetName();

    
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();
    const G4double kineticEnergy =
        step->GetPreStepPoint()->GetKineticEnergy();

    // Accumulate actual ionisation loss in the active GEM gas.  This is kept
    // separate from the truth boundary crossing below.
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.0) {
        if (volName == "DetectorIn1") fEventAction->AddEnergyIn1(edep);
        else if (volName == "DetectorIn2") fEventAction->AddEnergyIn2(edep);
        else if (volName == "DetectorOut1") fEventAction->AddEnergyOut1(edep);
        else if (volName == "DetectorOut2") fEventAction->AddEnergyOut2(edep);
        else if (volName == "DetectorOut3") fEventAction->AddEnergyOut3(edep);
        else if (volName == "DetectorTop") fEventAction->AddEnergyTop(edep);
    }

    // Truth positions and momenta are recorded only at volume entry.
    if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) return;
    const G4double time = step->GetPreStepPoint()->GetGlobalTime();

    if (volName == "DetectorIn1") fEventAction->SetHitIn1(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorIn2") fEventAction->SetHitIn2(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorOut0") fEventAction->SetHitOut0(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorOut1") fEventAction->SetHitOut1(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorOut2") fEventAction->SetHitOut2(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorOut3") fEventAction->SetHitOut3(pos, mom, time, kineticEnergy);
    else if (volName == "DetectorTop") fEventAction->SetHitTop(pos, mom, time, kineticEnergy);
    else if (volName == "PhysMound") {
        fEventAction->AddHitMound();
    }
    else if (volName == "PhysRoom") {
        fEventAction->AddHitRoom();
        fEventAction->AddHitRoomInterior();
    }
    else if (volName == "PhysRockWall") {
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
