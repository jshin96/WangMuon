#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void   EndOfRunAction(const G4Run* run);

    // Called once for every requested beam event, including events rejected
    // before Geant4 transport.  ROOT stores only four-layer coincidences, so
    // these counters are the authoritative generation-efficiency diagnostic.
    void RecordPrimaryGeneration(G4int eventID, G4long trials, G4bool accepted);

  private:
    G4long fGeneratedEvents = 0;
    G4long fAcceptedPrimaries = 0;
    G4long fAbortedPrimaries = 0;
    G4long fTotalTrials = 0;
    G4long fMaximumTrials = 0;
};

#endif
