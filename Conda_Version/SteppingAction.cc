#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const auto* particle = step->GetTrack()->GetDefinition();
    if (particle != G4MuonMinus::MuonMinusDefinition()
        && particle != G4MuonPlus::MuonPlusDefinition()) return;
    if (step->GetTrack()->GetParentID() != 0) return;
    if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) return;

    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;

    const G4String volName = volume->GetName();
    const auto position = step->GetPreStepPoint()->GetPosition();
    const auto momentum = step->GetPreStepPoint()->GetMomentum();
    const G4double radialDot = position.x() * momentum.x()
                             + position.y() * momentum.y();

    if (volName == "DetectorOuter") {
        if (radialDot < 0.0) fEventAction->SetHitIn1(position, momentum);
        else                 fEventAction->SetHitOut2(position, momentum);
    } else if (volName == "DetectorInner") {
        if (radialDot < 0.0) fEventAction->SetHitIn2(position, momentum);
        else                 fEventAction->SetHitOut1(position, momentum);
    }
}
