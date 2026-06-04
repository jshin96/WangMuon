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

    // Setters for the 4 layers
    void SetHitIn1 (G4ThreeVector pos, G4ThreeVector mom) { fPosIn1 = pos; fMomIn1 = mom; fHitIn1 = true; }
    void SetHitIn2 (G4ThreeVector pos, G4ThreeVector mom) { fPosIn2 = pos; fMomIn2 = mom; fHitIn2 = true; }
    void SetHitOut1(G4ThreeVector pos, G4ThreeVector mom) { fPosOut1 = pos; fMomOut1 = mom; fHitOut1 = true; }
    void SetHitOut2(G4ThreeVector pos, G4ThreeVector mom) { fPosOut2 = pos; fMomOut2 = mom; fHitOut2 = true; }

  private:
    G4ThreeVector fPosIn1, fMomIn1;
    G4ThreeVector fPosIn2, fMomIn2;
    G4ThreeVector fPosOut1, fMomOut1;
    G4ThreeVector fPosOut2, fMomOut2;
    
    G4bool fHitIn1, fHitIn2, fHitOut1, fHitOut2;
};

#endif
