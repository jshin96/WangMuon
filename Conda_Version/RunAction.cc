#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"

#include <cstdlib>

RunAction::RunAction() : G4UserRunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root"); // Or "csv"
    
    // Create the Ntuple
    analysisManager->CreateNtuple("MuonHits", "4-Layer Tracking Data");
    
    analysisManager->CreateNtupleDColumn("Trials");
    analysisManager->CreateNtupleIColumn("EventID");
    
    // Detector In 1 (Furthest out, -Y)
    analysisManager->CreateNtupleDColumn("In1_X");
    analysisManager->CreateNtupleDColumn("In1_Y");
    analysisManager->CreateNtupleDColumn("In1_Z");
    analysisManager->CreateNtupleDColumn("In1_Px");
    analysisManager->CreateNtupleDColumn("In1_Py");
    analysisManager->CreateNtupleDColumn("In1_Pz");

    // Detector In 2 (Closer to mound, -Y)
    analysisManager->CreateNtupleDColumn("In2_X");
    analysisManager->CreateNtupleDColumn("In2_Y");
    analysisManager->CreateNtupleDColumn("In2_Z");
    analysisManager->CreateNtupleDColumn("In2_Px");
    analysisManager->CreateNtupleDColumn("In2_Py");
    analysisManager->CreateNtupleDColumn("In2_Pz");

    // Detector Out 1 (Closer to mound, +Y)
    analysisManager->CreateNtupleDColumn("Out1_X");
    analysisManager->CreateNtupleDColumn("Out1_Y");
    analysisManager->CreateNtupleDColumn("Out1_Z");
    analysisManager->CreateNtupleDColumn("Out1_Px");
    analysisManager->CreateNtupleDColumn("Out1_Py");
    analysisManager->CreateNtupleDColumn("Out1_Pz");

    // Detector Out 2 (Furthest out, +Y)
    analysisManager->CreateNtupleDColumn("Out2_X");
    analysisManager->CreateNtupleDColumn("Out2_Y");
    analysisManager->CreateNtupleDColumn("Out2_Z");
    analysisManager->CreateNtupleDColumn("Out2_Px");
    analysisManager->CreateNtupleDColumn("Out2_Py");
    analysisManager->CreateNtupleDColumn("Out2_Pz");

    analysisManager->FinishNtuple();
}

RunAction::~RunAction() {
    delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*) {
    fGeneratedEvents = 0;
    fAcceptedPrimaries = 0;
    fAbortedPrimaries = 0;
    fTotalTrials = 0;
    fMaximumTrials = 0;
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("MoundTomographyData");
}

void RunAction::EndOfRunAction(const G4Run*) {
    const double meanTrials = fGeneratedEvents > 0
        ? static_cast<double>(fTotalTrials) / fGeneratedEvents : 0.0;
    G4cout << "PRIMARY_GENERATION_SUMMARY requested=" << fGeneratedEvents
           << " accepted=" << fAcceptedPrimaries
           << " aborted=" << fAbortedPrimaries
           << " mean_trials=" << meanTrials
           << " max_trials=" << fMaximumTrials << G4endl;
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

void RunAction::RecordPrimaryGeneration(G4int eventID, G4long trials,
                                        G4bool accepted) {
    ++fGeneratedEvents;
    fTotalTrials += trials;
    if (trials > fMaximumTrials) fMaximumTrials = trials;
    if (accepted) {
        ++fAcceptedPrimaries;
    } else {
        ++fAbortedPrimaries;
    }

    if (std::getenv("MOUND_PRINT_PRIMARY_DIRECTIONS") != nullptr) {
        G4cout << "PRIMARY_GENERATION event=" << eventID
               << " trials=" << trials
               << " status=" << (accepted ? "accepted" : "aborted")
               << G4endl;
    }
}
