#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"


void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());
    
    // Create the EventAction (Memory Bank)
    EventAction* eventAction = new EventAction();
    SetUserAction(eventAction);
    
    // Hand the memory bank to the SteppingAction
    SetUserAction(new SteppingAction(eventAction));
}
