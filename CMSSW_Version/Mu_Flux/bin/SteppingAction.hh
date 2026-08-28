#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction; // We only need to know that the note-taking class exists.

class SteppingAction : public G4UserSteppingAction {
public:
    // Give this walker a place where it can write detector notes.
    SteppingAction(EventAction* eventAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    EventAction* fEventAction; // The notebook for the current muon.
};

#endif
