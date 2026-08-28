#include "TrackingAction.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction(RunAction* runAction) 
: fRunAction(runAction) {}

void TrackingAction::PreUserTrackingAction(const G4Track*) {
    // We are currently handling all the muon data extraction 
    // inside SteppingAction.cc, so we can leave this blank!
}
