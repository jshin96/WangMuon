#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);

    // Setters called by SteppingAction
    void SetIncomingHit(G4ThreeVector pos, G4ThreeVector mom) {
        fPosIn = pos; 
        fMomIn = mom; 
        fHitIn = true;
    }
    
    void SetOutgoingHit(G4ThreeVector pos, G4ThreeVector mom) {
        fPosOut = pos; 
        fMomOut = mom; 
        fHitOut = true;
    }

  private:
    G4ThreeVector fPosIn, fMomIn;
    G4ThreeVector fPosOut, fMomOut;
    
    // Flags to ensure we only record muons that made it all the way through
    G4bool fHitIn;
    G4bool fHitOut;
};

#endif
