#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

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


    // -----------------------------------------------------------
    // FILTER 2: Did it just cross a boundary?
    // -----------------------------------------------------------
    // fGeomBoundary means the particle just entered a new volume.
    // If it is simply taking a step *inside* the same volume, ignore it.
    if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) {
        return; 
    }

    // -----------------------------------------------------------
    // FILTER 3: Which volume did it enter?
    // -----------------------------------------------------------
    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;

    G4String volName = volume->GetName();

    
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();

    // -----------------------------------------------------------
    // Transform the hit into the detector's own tilted-cylinder frame. A 2D
    // radial dot product there determines inward versus outward travel.
    // -----------------------------------------------------------
    const auto& globalToLocal =
        step->GetPreStepPoint()->GetTouchableHandle()
            ->GetHistory()->GetTopTransform();
    const G4ThreeVector localPos = globalToLocal.TransformPoint(pos);
    const G4ThreeVector localMom = globalToLocal.TransformAxis(mom);
    const double radialDot =
        localPos.x() * localMom.x() + localPos.y() * localMom.y();

    if (volName == "DetectorOuter") {
        if (radialDot < 0.0) {
            // Heading towards center = Initial Entry Plane
            fEventAction->SetHitIn1(pos, mom); 
        } else {
            // Heading away from center = Final Exit Plane
            fEventAction->SetHitOut2(pos, mom); 
        }
    } 
    else if (volName == "DetectorInner") {
        if (radialDot < 0.0) {
            // Heading towards center = Secondary Entry Plane
            fEventAction->SetHitIn2(pos, mom); 
        } else {
            // Heading away from center = Initial Exit Plane
            fEventAction->SetHitOut1(pos, mom); 
        }
    }
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
