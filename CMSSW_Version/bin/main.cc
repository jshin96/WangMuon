#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh" // Geant4's pre-packaged physics list

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv) {
    // 1. Create the Run Manager (The Physics Engine)
    G4RunManager* runManager = new G4RunManager;
    
    // 2. Register the THREE PILLARS
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT());
    runManager->SetUserInitialization(new ActionInitialization());

    runManager->Initialize();
    
    // 3. Get the pointer to the User Interface manager (To read the macro)
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    // 4. Force Batch Mode (Condor only)
    if (argc > 1) {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        G4cout << "Please provide a macro file (e.g., wang_muon run.mac)" << G4endl;
    }

    // 5. Clean up the engine
    delete runManager;
    return 0;
}
