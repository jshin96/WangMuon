#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
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
    
    // Define the boundaries of your sky (a 20m x 20m plane)
    fEcoMug->SetHSphereCenterPosition({{0., 0., 0}}); 
    fEcoMug->SetHSphereRadius(74.99 * m);

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
    // 1. Tell EcoMug to generate a new cosmic ray based on empirical data
    fEcoMug->Generate();
    double p_tot = fEcoMug->GetGenerationMomentum(); // in GeV

    // 2. Extract the physical properties from EcoMug
    // Note: EcoMug uses standard spherical angles naturally!
    /*
    std::array<double, 3> pos = fEcoMug->GetGenerationPosition();
    double theta = fEcoMug->GetGenerationTheta();
    double phi   = fEcoMug->GetGenerationPhi();
    // 3. Hand the EcoMug data over to the Geant4 Particle Gun
    fParticleGun->SetParticlePosition(G4ThreeVector(pos[0], pos[1], pos[2]));
    
    // Convert spherical angles to Cartesian momentum direction
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
        -std::abs(std::cos(theta)) // Ensure it always points downwards
    ));
    */

    // Testing with only local +Y direction muons. Rotate the source position
    // and direction with the rigid detector cylinder so the ray remains
    // parallel to the inclined ground and can cross all four detector layers.
    const char* inclineSetting = std::getenv("MOUND_GROUND_INCLINE_DEG");
    const G4double inclineDeg =
        inclineSetting ? std::atof(inclineSetting) : 15.0;
    const G4double incline = inclineDeg * CLHEP::pi / 180.0;
    const G4double inclineCos = std::cos(incline);
    const G4double inclineSin = std::sin(incline);
    const G4double localX = ((G4UniformRand()*15.0)-7.5)*m;
    const G4double localY = -74.95*m;
    const G4double axialHeight = G4UniformRand()*13.0*m;
    G4ThreeVector pos(localX,
                      inclineCos*localY + inclineSin*axialHeight,
                      -inclineSin*localY + inclineCos*axialHeight);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(
        G4ThreeVector(0., inclineCos, -inclineSin));
    
    // Set the kinetic energy
    // (For highly relativistic muons, p_tot ~ Kinetic Energy)
    fParticleGun->SetParticleEnergy(p_tot * GeV);

    // 4. Fire!
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
