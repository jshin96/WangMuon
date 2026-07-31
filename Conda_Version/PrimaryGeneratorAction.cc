#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ios.hh"
#include <array>
#include <cmath>
#include <cstdlib>

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    fParticleGun = new G4ParticleGun(1);

    // Set the default particle to a negative muon
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));

    // --- Initialize and Configure EcoMug ---
    fEcoMug = new EcoMug();
    
    // We want a standard half-spherical surface 
    fEcoMug->SetUseHSphere(); 
    
    // Keep the source hemisphere inside the 50 m world half-height while
    // enclosing the roughly 31 m detector/mound envelope.
    fEcoMug->SetHSphereCenterPosition({{0., 0., 0}}); 
    fEcoMug->SetHSphereRadius(45.0 * m);

    // Lock the momentum to your specific high-energy requirements (1 MeV to 100 GeV)
    // EcoMug calculates momentum in GeV/c by default
    fEcoMug->SetMinimumMomentum(10.);
    fEcoMug->SetMaximumMomentum(100.);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fEcoMug;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // EcoMug returns a point on the upper hemisphere and an inward cosmic
    // direction. Its returned theta already uses the downward convention.
    fEcoMug->Generate();
    const G4double momentumGeV = fEcoMug->GetGenerationMomentum();
    const std::array<double, 3>& position =
        fEcoMug->GetGenerationPosition();
    const G4double theta = fEcoMug->GetGenerationTheta();
    const G4double phi = fEcoMug->GetGenerationPhi();
    const G4ThreeVector direction(std::sin(theta)*std::cos(phi),
                                  std::sin(theta)*std::sin(phi),
                                  std::cos(theta));

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(
        particleTable->FindParticle(fEcoMug->GetCharge() < 0 ? "mu-" : "mu+"));
    fParticleGun->SetParticlePosition(
        G4ThreeVector(position[0], position[1], position[2]));
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticleMomentum(momentumGeV * GeV);

    if (std::getenv("MOUND_PRINT_PRIMARY_DIRECTIONS") != nullptr) {
        G4cout << "COSMIC_PRIMARY event=" << anEvent->GetEventID()
               << " position_m=(" << position[0]/m << ","
               << position[1]/m << "," << position[2]/m << ")"
               << " direction=(" << direction.x() << ","
               << direction.y() << "," << direction.z() << ")"
               << " momentum_GeV=" << momentumGeV
               << " charge=" << fEcoMug->GetCharge() << G4endl;
    }

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
