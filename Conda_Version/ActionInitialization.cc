#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"      
#include "EventAction.hh"     
#include "SteppingAction.hh"  

ActionInitialization::ActionInitialization()
 : G4VUserActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

// This is called in multithreading mode to merge the ROOT files
void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

// This is where the worker threads get their action classes
void ActionInitialization::Build() const
{
    // 1. Run Action owns generation-efficiency counters.
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);

    // 2. Particle Gun reports each requested primary to RunAction.
    SetUserAction(new PrimaryGeneratorAction(runAction));

    // 3. Event Action (handles the coincidence trigger)
    EventAction* eventAction = new EventAction();
    SetUserAction(eventAction);

    // 4. Stepping Action (needs the EventAction pointer to pass the hits)
    SetUserAction(new SteppingAction(eventAction));
}
