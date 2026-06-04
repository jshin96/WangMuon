#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // Get the volume of the current step
    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!volume) return;

    G4String volName = volume->GetName();

    // Route the data based on which detector was hit
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
}
