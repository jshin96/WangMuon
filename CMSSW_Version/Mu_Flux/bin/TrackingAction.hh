#ifndef TRACKING_ACTION_HH
#define TRACKING_ACTION_HH

#include "G4UserTrackingAction.hh"
#include "RunAction.hh"

// A spare hook Geant4 offers when a new particle track starts.
class TrackingAction : public G4UserTrackingAction {
public:
    TrackingAction(RunAction* runAction);
    ~TrackingAction() override = default;
    void PreUserTrackingAction(const G4Track* track) override;

private:
    RunAction* fRunAction;
};
#endif
