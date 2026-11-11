#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);

    // Truth is the boundary crossing; GEM data are produced independently at
    // end of event from gas ionisation, gain, electronics noise and threshold.
    void SetHitIn1 (G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitIn2 (G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitOut0(G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitOut1(G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitOut2(G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitOut3(G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void SetHitTop (G4ThreeVector pos, G4ThreeVector mom, G4double time,
                    G4double kineticEnergy);
    void AddEnergyIn1(G4double energy);
    void AddEnergyIn2(G4double energy);
    void AddEnergyOut1(G4double energy);
    void AddEnergyOut2(G4double energy);
    void AddEnergyOut3(G4double energy);
    void AddEnergyTop(G4double energy);

    void AddHitMound() { fHitMound = true; }
    void AddHitRoom()  { fHitRoom = true; }
    void AddHitRoomInterior() { fHitRoomInterior = true; }
  private:
  public: // Used only by the ROOT-output helper in EventAction.cc.
    struct Hit {
        G4ThreeVector pos, mom, gemPos;
        G4double time, kineticEnergy, energy, gemTime, gemCharge;
        G4bool hit, gemValid;
    };
  private:
    void SetTruth(Hit& hit, G4ThreeVector pos, G4ThreeVector mom, G4double time,
                  G4double kineticEnergy);
    void Reset(Hit& hit);
    void Digitize(Hit& hit, const G4ThreeVector& firstAxis,
                  const G4ThreeVector& secondAxis);
    Hit fIn1, fIn2, fOut0, fOut1, fOut2, fOut3, fTop;
    G4ThreeVector fGEMTangentialAxis, fGEMAxialAxis;
    G4ThreeVector fTopFirstAxis, fTopSecondAxis;
    RunAction* fRunAction;
    bool fNeutronMode;
    bool fTopNormalizationMode;
    bool fHitMound;
    bool fHitRoom;
    bool fHitRoomInterior;
};

#endif
