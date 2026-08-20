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

    // Truth is the boundary crossing; GEM data are produced independently at
    // end of event from gas ionisation, gain, electronics noise and threshold.
    void SetHitIn1 (G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void SetHitIn2 (G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void SetHitOut1(G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void SetHitOut2(G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void SetHitOut3(G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void AddEnergyIn1(G4double energy);
    void AddEnergyIn2(G4double energy);
    void AddEnergyOut1(G4double energy);
    void AddEnergyOut2(G4double energy);
    void AddEnergyOut3(G4double energy);

    void AddHitMound() { fHitMound = true; }
    void AddHitRoom()  { fHitRoom = true; }
  private:
  public: // Used only by the ROOT-output helper in EventAction.cc.
    struct Hit {
        G4ThreeVector pos, mom, gemPos;
        G4double time, energy, gemTime, gemCharge;
        G4bool hit, gemValid;
    };
  private:
    void SetTruth(Hit& hit, G4ThreeVector pos, G4ThreeVector mom, G4double time);
    void Reset(Hit& hit);
    void Digitize(Hit& hit);
    Hit fIn1, fIn2, fOut1, fOut2, fOut3;
    G4ThreeVector fGEMTangentialAxis, fGEMAxialAxis;
    bool fHitMound;
    bool fHitRoom;
};

#endif
