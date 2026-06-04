#include "RunAction.hh"
#include "G4AnalysisManager.hh"

RunAction::RunAction() : G4UserRunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root"); // Or "csv"
    
    // Create the Ntuple
    analysisManager->OpenFile("MoundTomographyData");
    analysisManager->CreateNtuple("MuonHits", "4-Layer Tracking Data");
    
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
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}
