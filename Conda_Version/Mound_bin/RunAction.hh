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
    // before Geant4 transport. ROOT stores only IN2-triggered rows in side
    // mode and top-monitor-triggered rows in normalization mode, so these
    // counters are the authoritative generation-efficiency diagnostic.
    void RecordPrimaryGeneration(G4int eventID, G4long trials, G4bool accepted);
    void SetSourceFluxMetadata(G4double ratePerAreaHzM2,
                               G4double ratePerAreaErrorHzM2,
                               G4double sourceAreaM2,
                               G4long integrationPoints);
    // 0 = side-cylinder cosmic muons, 1 = top-sky cosmic muons,
    // 2 = directed flat-energy neutron beam.
    void SetSourceMode(G4int sourceMode);
    void SetNeutronEnergyRange(G4double minimumGeV, G4double maximumGeV);
    void RecordFiveGEMEvent();
    void RecordTopTruthEvent();
    void RecordTopGEMEvent();
    void RecordTopRoomEvent();

  private:
    G4long fGeneratedEvents = 0;
    G4long fAcceptedPrimaries = 0;
    G4long fAbortedPrimaries = 0;
    G4long fTotalTrials = 0;
    G4long fMaximumTrials = 0;
    G4long fFiveGEMEvents = 0;
    G4long fTopTruthEvents = 0;
    G4long fTopGEMEvents = 0;
    G4long fTopRoomEvents = 0;
    G4long fRateIntegrationPoints = 0;
    G4double fSourceRatePerAreaHzM2 = 0.0;
    G4double fSourceRatePerAreaErrorHzM2 = 0.0;
    G4double fSourceAreaM2 = 0.0;
    G4int fSourceMode = 0;
    G4double fNeutronMinKineticEnergyGeV = 0.0;
    G4double fNeutronMaxKineticEnergyGeV = 0.0;
};

#endif
