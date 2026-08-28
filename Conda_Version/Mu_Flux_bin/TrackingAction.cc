#include "TrackingAction.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"

TrackingAction::TrackingAction(RunAction* runAction) 
: fRunAction(runAction) {}

void TrackingAction::PreUserTrackingAction(const G4Track*) {
    // Geant4 discards trajectories unless explicitly asked to retain them.
    // The visualization scene can then render each retained trajectory after
    // a single /run/beamOn event.
    fpTrackingManager->SetStoreTrajectory(1);
}
