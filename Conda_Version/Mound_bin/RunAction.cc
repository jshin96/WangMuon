#include <cmath>
#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"

RunAction::RunAction() : G4UserRunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root"); // Or "csv"
    analysisManager->CreateNtuple("MuonHits", "Side/top truth and parameterised GEM data");
    
    // PrimaryGeneratorAction records the number of EcoMug attempts with
    // FillNtupleDColumn, so this must be a double column.  Keeping the types
    // matched prevents one Geant4 warning (and log write) per accepted event.
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
    // Zenith relative to the downward vertical (-Z), evaluated from the
    // incoming truth momentum at IN2.  Units are degrees.
    analysisManager->CreateNtupleDColumn("IncomingZenithAngle_deg");

    // Air-only truth-exit slice just beyond the mound.  It deliberately has
    // no momentum, GEM, or trigger-validity fields.
    analysisManager->CreateNtupleDColumn("Out0_X");
    analysisManager->CreateNtupleDColumn("Out0_Y");
    analysisManager->CreateNtupleDColumn("Out0_Z");

    // Detector Out 1 (Closer to mound, +Y)
    analysisManager->CreateNtupleDColumn("Out1_X");
    analysisManager->CreateNtupleDColumn("Out1_Y");
    analysisManager->CreateNtupleDColumn("Out1_Z");
    analysisManager->CreateNtupleDColumn("Out1_Px");
    analysisManager->CreateNtupleDColumn("Out1_Py");
    analysisManager->CreateNtupleDColumn("Out1_Pz");

    // Detector Out 2 (second pre-tungsten +Y layer)
    analysisManager->CreateNtupleDColumn("Out2_X");
    analysisManager->CreateNtupleDColumn("Out2_Y");
    analysisManager->CreateNtupleDColumn("Out2_Z");
    analysisManager->CreateNtupleDColumn("Out2_Px");
    analysisManager->CreateNtupleDColumn("Out2_Py");
    analysisManager->CreateNtupleDColumn("Out2_Pz");

    // Detector Out 3 (post-tungsten +Y layer)
    analysisManager->CreateNtupleDColumn("Out3_X");
    analysisManager->CreateNtupleDColumn("Out3_Y");
    analysisManager->CreateNtupleDColumn("Out3_Z");
    analysisManager->CreateNtupleDColumn("Out3_Px");
    analysisManager->CreateNtupleDColumn("Out3_Py");
    analysisManager->CreateNtupleDColumn("Out3_Pz");

    analysisManager->CreateNtupleIColumn("TrajectoryFlag");

    analysisManager->CreateNtupleIColumn("Out1Valid");
    analysisManager->CreateNtupleIColumn("Out2Valid");
    analysisManager->CreateNtupleIColumn("Out3Valid");

    // GEM digitisation.  Coordinates/times are reconstructed measurements;
    // Edep is Geant4 gas truth and Charge contains gain/noise/threshold effects.
    const char* planes[] = {"GEMIn1", "GEMIn2", "GEMOut1", "GEMOut2", "GEMOut3"};
    for (const char* plane : planes) {
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_X");
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_Y");
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_Z");
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_Time_ns");
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_Edep_keV");
        analysisManager->CreateNtupleDColumn(G4String(plane) + "_Charge_fC");
        analysisManager->CreateNtupleIColumn(G4String(plane) + "_Valid");
    }
    // Optional horizontal monitor at the mound top. These branches remain
    // present (and NaN/invalid) in the default side-cylinder mode so both
    // source modes retain one stable ROOT schema.
    analysisManager->CreateNtupleDColumn("Top_X");
    analysisManager->CreateNtupleDColumn("Top_Y");
    analysisManager->CreateNtupleDColumn("Top_Z");
    analysisManager->CreateNtupleDColumn("Top_Px");
    analysisManager->CreateNtupleDColumn("Top_Py");
    analysisManager->CreateNtupleDColumn("Top_Pz");
    analysisManager->CreateNtupleDColumn("TopIncomingZenithAngle_deg");
    analysisManager->CreateNtupleIColumn("TopValid");
    analysisManager->CreateNtupleIColumn("ReachedRoomInterior");
    analysisManager->CreateNtupleDColumn("GEMTop_X");
    analysisManager->CreateNtupleDColumn("GEMTop_Y");
    analysisManager->CreateNtupleDColumn("GEMTop_Z");
    analysisManager->CreateNtupleDColumn("GEMTop_Time_ns");
    analysisManager->CreateNtupleDColumn("GEMTop_Edep_keV");
    analysisManager->CreateNtupleDColumn("GEMTop_Charge_fC");
    analysisManager->CreateNtupleIColumn("GEMTop_Valid");
    // Additive primary-particle fields. Keeping these at the end preserves
    // every pre-existing column index used by the cosmic-muon analysis.
    analysisManager->CreateNtupleIColumn("PrimaryPDG");
    analysisManager->CreateNtupleDColumn("In2_KineticEnergy_GeV");
    analysisManager->CreateNtupleDColumn("Out1_KineticEnergy_GeV");
    analysisManager->FinishNtuple();

    // One row per run.  Source quantities refer to the fully configured
    // EcoMug phase space; acceptance rates are normalized by all trials.
    analysisManager->CreateNtuple("RunMetadata", "Flux normalization and exposure summary");
    analysisManager->CreateNtupleDColumn("SourceRatePerArea_Hz_m2");
    analysisManager->CreateNtupleDColumn("SourceRatePerAreaError_Hz_m2");
    analysisManager->CreateNtupleDColumn("SourceArea_m2");
    analysisManager->CreateNtupleDColumn("SourceRate_Hz");
    analysisManager->CreateNtupleIColumn("RateIntegrationPoints");
    analysisManager->CreateNtupleIColumn("RequestedEvents");
    analysisManager->CreateNtupleIColumn("AcceptedIn1In2");
    analysisManager->CreateNtupleIColumn("AbortedPrimaries");
    analysisManager->CreateNtupleDColumn("TotalTrials");
    analysisManager->CreateNtupleIColumn("RecordedFiveGEM");
    analysisManager->CreateNtupleDColumn("In1In2Rate_Hz");
    analysisManager->CreateNtupleDColumn("RecordedFiveGEMRate_Hz");
    analysisManager->CreateNtupleDColumn("EquivalentLiveTime_s");
    analysisManager->CreateNtupleDColumn("DaysFor10000FiveGEM");
    analysisManager->CreateNtupleIColumn("SourceMode"); // 0: side muon, 1: top muon, 2: neutron
    analysisManager->CreateNtupleIColumn("RecordedTopTruth");
    analysisManager->CreateNtupleIColumn("RecordedTopGEM");
    analysisManager->CreateNtupleIColumn("TopReachedRoomInterior");
    analysisManager->CreateNtupleDColumn("TopTruthRate_Hz");
    analysisManager->CreateNtupleDColumn("TopGEMRate_Hz");
    analysisManager->CreateNtupleDColumn("TopRoomRate_Hz");
    analysisManager->CreateNtupleDColumn("DaysFor10000TopGEM");
    analysisManager->CreateNtupleDColumn("TopRoomTransmission");
    analysisManager->CreateNtupleDColumn("TopRoomTransmissionError");
    analysisManager->CreateNtupleDColumn("NeutronMinKineticEnergy_GeV");
    analysisManager->CreateNtupleDColumn("NeutronMaxKineticEnergy_GeV");
    analysisManager->FinishNtuple();
}

RunAction::~RunAction() {
}

void RunAction::BeginOfRunAction(const G4Run*) {
    fGeneratedEvents = 0;
    fAcceptedPrimaries = 0;
    fAbortedPrimaries = 0;
    fTotalTrials = 0;
    fMaximumTrials = 0;
    fFiveGEMEvents = 0;
    fTopTruthEvents = 0;
    fTopGEMEvents = 0;
    fTopRoomEvents = 0;
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
    const G4double sourceRateHz = fSourceRatePerAreaHzM2*fSourceAreaM2;
    const G4double sourceRateErrorHz = fSourceRatePerAreaErrorHzM2*fSourceAreaM2;
    const G4double trialFraction = fTotalTrials > 0
        ? 1.0/static_cast<G4double>(fTotalTrials) : 0.0;
    const G4double in1In2RateHz = sourceRateHz*fAcceptedPrimaries*trialFraction;
    const G4double recordedRateHz = sourceRateHz*fFiveGEMEvents*trialFraction;
    const G4double topTruthRateHz = sourceRateHz*fTopTruthEvents*trialFraction;
    const G4double topGEMRateHz = sourceRateHz*fTopGEMEvents*trialFraction;
    const G4double topRoomRateHz = sourceRateHz*fTopRoomEvents*trialFraction;
    const G4double equivalentLiveTime = sourceRateHz > 0.0
        ? static_cast<G4double>(fTotalTrials)/sourceRateHz : 0.0;
    const G4double daysFor10000 = recordedRateHz > 0.0
        ? 10000.0/(recordedRateHz*86400.0) : 0.0;
    const G4double daysFor10000Top = topGEMRateHz > 0.0
        ? 10000.0/(topGEMRateHz*86400.0) : 0.0;
    const G4double topRoomTransmission = fTopTruthEvents > 0
        ? static_cast<G4double>(fTopRoomEvents)/fTopTruthEvents : 0.0;
    const G4double topRoomTransmissionError = fTopTruthEvents > 0
        ? std::sqrt(topRoomTransmission*(1.0 - topRoomTransmission)
                    /static_cast<G4double>(fTopTruthEvents))
        : 0.0;

    constexpr G4int metadataNtuple = 1;
    G4int column = 0;
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, fSourceRatePerAreaHzM2);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, fSourceRatePerAreaErrorHzM2);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, fSourceAreaM2);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, sourceRateHz);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fRateIntegrationPoints);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fGeneratedEvents);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fAcceptedPrimaries);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fAbortedPrimaries);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, static_cast<G4double>(fTotalTrials));
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fFiveGEMEvents);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, in1In2RateHz);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, recordedRateHz);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, equivalentLiveTime);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, daysFor10000);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fSourceMode);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fTopTruthEvents);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fTopGEMEvents);
    analysisManager->FillNtupleIColumn(metadataNtuple, column++, fTopRoomEvents);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, topTruthRateHz);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, topGEMRateHz);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, topRoomRateHz);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, daysFor10000Top);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, topRoomTransmission);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++, topRoomTransmissionError);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++,
                                       fNeutronMinKineticEnergyGeV);
    analysisManager->FillNtupleDColumn(metadataNtuple, column++,
                                       fNeutronMaxKineticEnergyGeV);
    analysisManager->AddNtupleRow(metadataNtuple);

    G4cout << "RUN_METADATA source_rate_hz=" << sourceRateHz
           << " +/- " << sourceRateErrorHz
           << " total_trials=" << fTotalTrials
           << " recorded_five_gem=" << fFiveGEMEvents
           << " recorded_top_gem=" << fTopGEMEvents
           << " top_reached_room=" << fTopRoomEvents
           << " top_room_transmission=" << topRoomTransmission
           << " +/- " << topRoomTransmissionError
           << " recorded_rate_per_day=" << recordedRateHz*86400.0
           << " equivalent_live_time_s=" << equivalentLiveTime
           << " days_for_10000_five_gem=" << daysFor10000 << G4endl;
    analysisManager->Write();
    analysisManager->CloseFile();
}

void RunAction::SetSourceFluxMetadata(G4double ratePerAreaHzM2,
                                      G4double ratePerAreaErrorHzM2,
                                      G4double sourceAreaM2,
                                      G4long integrationPoints) {
    fSourceRatePerAreaHzM2 = ratePerAreaHzM2;
    fSourceRatePerAreaErrorHzM2 = ratePerAreaErrorHzM2;
    fSourceAreaM2 = sourceAreaM2;
    fRateIntegrationPoints = integrationPoints;
}

void RunAction::SetSourceMode(G4int sourceMode) {
    fSourceMode = sourceMode;
}

void RunAction::SetNeutronEnergyRange(G4double minimumGeV,
                                      G4double maximumGeV) {
    fNeutronMinKineticEnergyGeV = minimumGeV;
    fNeutronMaxKineticEnergyGeV = maximumGeV;
}

void RunAction::RecordFiveGEMEvent() {
    ++fFiveGEMEvents;
}

void RunAction::RecordTopTruthEvent() {
    ++fTopTruthEvents;
}

void RunAction::RecordTopGEMEvent() {
    ++fTopGEMEvents;
}

void RunAction::RecordTopRoomEvent() {
    ++fTopRoomEvents;
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

    G4cout << "PRIMARY_GENERATION event=" << eventID
           << " trials=" << trials
           << " status=" << (accepted ? "accepted" : "aborted")
           << G4endl;
}
