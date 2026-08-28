#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
// This class makes one muon for every event.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  // Put a new muon at the start of the beam line.
  void GeneratePrimaries(G4Event*) override;

private:
  G4ParticleGun* fGun;  // The particle-making tool.
  G4double fSourceX;
  G4double fHalfBeamSize;
  G4double fMinEnergy;
  G4double fMaxEnergy;
  G4double fBeamTheta;
  G4double fBeamPhi;
};
#endif
