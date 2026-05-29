#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include "G4UserSteppingAction.hh"
#include "RunAction.hh"

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(RunAction* runAction);
    ~SteppingAction() override = default;
    void UserSteppingAction(const G4Step* step) override;

private:
    RunAction* fRunAction;
};
#endif
