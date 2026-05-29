#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RunManager.hh"
#include "FTFP_BERT.hh" // Geant4's pre-packaged physics list

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv) {
    // 1. Create the User Interface
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    // 2. Create the Run Manager
//    auto* runManager = G4RunManagerFactory::CreateRunManager();
    G4RunManager* runManager = new G4RunManager;
    // 3. Register the THREE PILLARS
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT());
    runManager->SetUserInitialization(new ActionInitialization());

    runManager->Initialize();

    // 4. Initialize Visualization
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    // 5. Tell the UI to draw our box
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");

    // 6. Start the interactive session!
    ui->SessionStart();

    // Clean up after the window is closed
    delete ui;
    delete visManager;
    delete runManager;
    return 0;
}
