#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include "G4UserEventAction.hh"
#include "RunAction.hh"

class EventAction : public G4UserEventAction {
public:
    // We pass the RunAction pointer so this class can access the ClearVectors() method
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

private:
    RunAction* fRunAction;
};
#endif
