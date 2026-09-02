#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "FTFP_BERT_HP.hh" // Geant4's pre-packaged physics list
#include "G4RunManagerFactory.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {
void ConfigureRandomSeed() {
    const char* value = std::getenv("MUON_RANDOM_SEED");
    if (!value || *value == '\0') {
        std::cout << "Geant4 random seed: default engine state" << std::endl;
        return;
    }

    char* end = nullptr;
    errno = 0;
    const unsigned long seed = std::strtoul(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0' || seed == 0 ||
        seed > static_cast<unsigned long>(std::numeric_limits<long>::max())) {
        std::cerr << "Invalid MUON_RANDOM_SEED='" << value
                  << "'; expected a positive integer accepted by CLHEP."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    G4Random::setTheSeed(static_cast<long>(seed));
    std::cout << "Geant4 random seed: " << seed << std::endl;
}
}  // namespace



int main(int argc, char** argv) {
    // Condor supplies a different MUON_RANDOM_SEED per job.  Do this before
    // Geant4 initializes physics so all random draws belong to this stream.
    ConfigureRandomSeed();
    // The run manager is the conductor that tells all Geant4 helpers when to work.
    G4RunManager* runManager = new G4RunManager;



    // Tell Geant4 how to build the room, which physics rules to use, and what
    // the helpers should do when a muon is made or moves.
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new FTFP_BERT_HP());
    runManager->SetUserInitialization(new ActionInitialization());

    runManager->Initialize();

    // Read the run.mac instructions in batch jobs.
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // 4. Batch mode for production; interactive UI with a drawn geometry when
    // no macro is supplied.
    if (argc > 1) {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        // Never construct a visualisation manager in batch/Condor mode.
        // The available interactive driver is selected by this Geant4 build.
        auto* visManager = new G4VisExecutive;
        visManager->Initialize();
        auto* ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/vis/open OGL");
        UImanager->ApplyCommand("/vis/drawVolume");
        UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
        UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate 100");
        ui->SessionStart();
        delete ui;
        delete visManager;
    }


    // 5. Clean up the engine
    delete runManager;
    return 0;

} 
