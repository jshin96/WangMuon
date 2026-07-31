#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "EcoMug.h" // <-- Include the new library

class G4Event;
class RunAction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(RunAction* runAction);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun* fParticleGun;
    EcoMug* fEcoMug;
    RunAction* fRunAction;
};

#endif
