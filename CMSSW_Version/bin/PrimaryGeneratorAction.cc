#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    fParticleGun = new G4ParticleGun(1);

    // Set the default particle to a negative muon
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));

    // --- Initialize and Configure EcoMug ---
    fEcoMug = new EcoMug();
    
    // We want a standard flat surface generation (a plane above your dirt mound)
    fEcoMug->SetUseSky(); 
    
    // Define the boundaries of your sky (a 20m x 20m plane)
    fEcoMug->SetSkyCenterPosition({{0., 0., 20.1 * m}}); 
    fEcoMug->SetSkySize({{60. * m, 60. * m}});

    // Lock the momentum to your specific high-energy requirements (1 MeV to 100 GeV)
    // EcoMug calculates momentum in GeV/c by default
    fEcoMug->SetMinimumMomentum(0.001);
    fEcoMug->SetMaximumMomentum(100.);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fEcoMug;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // 1. Tell EcoMug to generate a new cosmic ray based on empirical data
    fEcoMug->Generate();

    // 2. Extract the physical properties from EcoMug
    // Note: EcoMug uses standard spherical angles naturally!
    std::array<double, 3> pos = fEcoMug->GetGenerationPosition();
    double p_tot = fEcoMug->GetGenerationMomentum(); // in GeV
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
    
    // Set the kinetic energy
    // (For highly relativistic muons, p_tot ~ Kinetic Energy)
    fParticleGun->SetParticleEnergy(p_tot * GeV);

    // 4. Fire!
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
