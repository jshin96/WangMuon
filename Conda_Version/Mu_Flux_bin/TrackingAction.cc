#include "TrackingAction.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction(RunAction* runAction) 
: fRunAction(runAction) {}

void TrackingAction::PreUserTrackingAction(const G4Track*) {
    // SteppingAction writes the detector notes, so this hook has no job yet.
}
