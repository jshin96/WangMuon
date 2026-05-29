#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction());
    
    // 1. Create the RunAction first
    auto runAction = new RunAction();
    SetUserAction(runAction);

    // 2. Pass its pointer to the others so they all share the exact same vectors
    SetUserAction(new EventAction(runAction));
    SetUserAction(new TrackingAction(runAction));
    SetUserAction(new SteppingAction(runAction));
}
