#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
	auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	if (volume->GetName() == "DetectorMinusY") {
	    fEventAction->SetIncomingHit(
		step->GetPreStepPoint()->GetPosition(),
		step->GetPreStepPoint()->GetMomentum()
	    );
	}
	else if (volume->GetName() == "DetectorPlusY") {
	    fEventAction->SetOutgoingHit(
		step->GetPreStepPoint()->GetPosition(),
		step->GetPreStepPoint()->GetMomentum()
	    );
	}
}
