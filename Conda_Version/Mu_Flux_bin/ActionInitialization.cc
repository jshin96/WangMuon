#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"      
#include "EventAction.hh"     
#include "SteppingAction.hh"  
#include "TrackingAction.hh"

ActionInitialization::ActionInitialization()
 : G4VUserActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

// The master only needs the file helper; workers make and follow muons.
void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

// Give each worker all of the jobs it needs for one event.
void ActionInitialization::Build() const
{
    RunAction* runAction = new RunAction(); // Makes the ROOT output file.
    SetUserAction(runAction);
    SetUserAction(new PrimaryGeneratorAction()); // Makes the +X beam muon.
    EventAction* eventAction = new EventAction(); // Holds its three GEM notes.
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(eventAction));
    SetUserAction(new TrackingAction(runAction));
}
