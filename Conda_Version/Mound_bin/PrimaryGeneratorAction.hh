#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "EcoMug.h" // <-- Include the new library

class G4Event;
class RunAction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(RunAction* runAction);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
    struct FlatPlane {
        G4double width;
        G4double height;
        G4double y;
        G4double z;
    };

    G4ParticleGun* fParticleGun;
    EcoMug* fEcoMug;
    RunAction* fRunAction;
    G4double fGroundSlope;
    G4double fMoundRadius;
    FlatPlane fIn1;
    FlatPlane fIn2;
    G4long fMaximumAttempts;
    G4bool fTopNormalizationMode;
    G4double fTopWidth;
    G4double fTopLength;
    G4ThreeVector fTopCenter;
};

#endif
