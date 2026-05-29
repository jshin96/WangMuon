#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"



PrimaryGeneratorAction::PrimaryGeneratorAction() {
    fParticleGun = new G4ParticleGun(1); // Shoot 1 particle per event

    // Define the particle as a negative muon
    auto particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
    
    // Shoot it straight down the Z-axis with 100 MeV of energy
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
//    fParticleGun->SetParticleEnergy(100.*MeV);
    
    // Start it at z = -40 cm (just outside our 10cm water box)
//    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -40.*cm));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // 1. Random Energy 
    double energy = (10. + G4UniformRand() * 90.) * GeV;
    fParticleGun->SetParticleEnergy(energy);

    // 2. Position: Uniform on a 20m disk above the mound
    double r = 20.*m;
    // Only initiate in any X, only negative Y with 1/4 pi overhang, and positive Z with sin distribution of probability
    double posTheta = std::asin(G4UniformRand());
//    double posPhi = CLHEP::pi * (1.5 * G4UniformRand() - 0.25) ;
//    fParticleGun->SetParticlePosition(G4ThreeVector(r * std::cos(posPhi) * std::sin(posTheta), -r * std::sin(posPhi) * std::sin(posTheta), r * std::cos(posTheta)));
    double posX = 30*(G4UniformRand()-0.5)  ;
    double posY = 50*m  ;
    fParticleGun->SetParticlePosition(G4ThreeVector(posX*m, posY, 1*m));

    // 3. Direction: Zenith angle weighted by cos(theta)
    // std::asin(G4UniformRand()) mathematically gives a cos(theta) weighted distribution
    double theta = std::asin(G4UniformRand()); 
    double dirPhi = CLHEP::pi * G4UniformRand();
    // Particle shooting down with Pz with sine weight, only in positive Y direction, 
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(
//        std::sin(theta) * std::cos(dirPhi),
//        std::sin(theta) * std::sin(dirPhi),
//        -std::cos(theta) // Negative Z points it straight down into the dirt
//    ));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(
        0,
        -1,
        0 
    ));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
