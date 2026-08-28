#include "RunAction.hh"
#include "G4AnalysisManager.hh"
RunAction::RunAction(){
  // Make a ROOT table.  Each row is one muon that reached all three GEMs.
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("root");
  analysis->CreateNtuple("MuonHits", "Three-GEM beam spectrometer");
  analysis->CreateNtupleIColumn("EventID");
  const char* planes[] = {"GEMIn1", "GEMIn2", "GEMOut"};
  for (const auto* plane : planes) {
    // Truth is Geant4's exact answer; reco is the fuzzy detector answer.
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthX_mm");
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthY_mm");
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthZ_mm");
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthPx_GeV");
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthPy_GeV");
    analysis->CreateNtupleDColumn(G4String(plane) + "_TruthPz_GeV");
    analysis->CreateNtupleDColumn(G4String(plane) + "_RecoX_mm");
    analysis->CreateNtupleDColumn(G4String(plane) + "_RecoY_mm");
    analysis->CreateNtupleDColumn(G4String(plane) + "_RecoZ_mm");
  }
  analysis->CreateNtupleDColumn("TruthDeflectionAngle_rad");
  analysis->CreateNtupleDColumn("RecoDeflectionChordAngle_rad");
  analysis->FinishNtuple();
}
void RunAction::BeginOfRunAction(const G4Run*){ // Open the file before the first muon.
  G4AnalysisManager::Instance()->OpenFile("MuonBeamSpectrometer");
}
void RunAction::EndOfRunAction(const G4Run*){ // Write the finished notebook to disk.
  auto* analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile();
}
