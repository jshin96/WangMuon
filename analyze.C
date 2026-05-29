#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

void analyze() {
    // 1. Open the file and grab the tree
    TFile *file = TFile::Open("build/wang_muon_data.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    
    TTree *tree = (TTree*)file->Get("ParticleTree");

    // 2. Define pointers to vectors and initialize them to nullptr
    std::vector<int> *particle_ID = nullptr;
    std::vector<int> *particle_Mother = nullptr;
    std::vector<int> *particle_PID = nullptr;
    std::vector<double> *particle_Energy = nullptr;
    
    // (You can do the same for the Step vectors here)
    std::vector<double> *step_PosX = nullptr;
    std::vector<double> *step_PosY = nullptr;
    std::vector<double> *step_PosZ = nullptr;

    // 3. Link the tree branches to our pointers
    tree->SetBranchAddress("Particle_ID", &particle_ID);
    tree->SetBranchAddress("Particle_Mother", &particle_Mother);
    tree->SetBranchAddress("Particle_PID", &particle_PID);
    tree->SetBranchAddress("Particle_Energy", &particle_Energy);
    
    tree->SetBranchAddress("Step_PosX", &step_PosX);
    tree->SetBranchAddress("Step_PosY", &step_PosY);
    tree->SetBranchAddress("Step_PosZ", &step_PosZ);

    // 4. Loop over every event (entry) in the tree
    int nEvents = tree->GetEntries();
    std::cout << "Total Events to process: " << nEvents << std::endl;

    for (int i = 0; i < nEvents; i++) {
        // Load the data for event 'i' into our vectors
        tree->GetEntry(i);
        
        std::cout << "\n=== Event " << i << " ===" << std::endl;
        std::cout << "Total particles generated: " << particle_ID->size() << std::endl;
        std::cout << "Total steps recorded: " << step_PosX->size() << std::endl;

        // 5. Loop over the particles INSIDE this specific event
        int numPhotons = 0;
        for (size_t j = 0; j < particle_ID->size(); j++) {
            
            // If Mother == 0, this is our primary particle from the gun!
            if (particle_Mother->at(j) == 0) {
                std::cout << "  -> Primary Particle: PID = " << particle_PID->at(j) 
                          << " | Initial Energy = " << particle_Energy->at(j) << " MeV" << std::endl;
            }
            
            // Count secondary photons (PDG ID for photon is 22)
            if (particle_PID->at(j) == 22) {
                numPhotons++;
            }
        }
        
        std::cout << "  -> Secondary photons created in this shower: " << numPhotons << std::endl;
    }

    // Clean up
    file->Close();
}
