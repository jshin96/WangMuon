#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction; // Forward declaration

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(EventAction* eventAction); // Require the memory bank
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    EventAction* fEventAction; // Store the pointer
};

#endif
