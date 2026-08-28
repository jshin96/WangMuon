#include "MoundTomographyEMML_MPI_TV.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH3D.h"
#include "TH1D.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <mpi.h>

namespace {
constexpr double kMoundRadiusMm = 26500.0;
constexpr double kMoundHeightMm = 12700.0;
bool IntersectSlopedMound(const TVector3& p, const TVector3& d, double slope,
                          bool forward, TVector3& intersection) {
    const double zr = p.Z() + slope*p.Y(), dzr = d.Z() + slope*d.Y();
    const double r2 = kMoundRadiusMm*kMoundRadiusMm, h2 = kMoundHeightMm*kMoundHeightMm;
    const double a = (d.X()*d.X()+d.Y()*d.Y())/r2 + dzr*dzr/h2;
    const double b = 2.0*((p.X()*d.X()+p.Y()*d.Y())/r2 + zr*dzr/h2);
    const double c = (p.X()*p.X()+p.Y()*p.Y())/r2 + zr*zr/h2 - 1.0;
    const double disc = b*b - 4.0*a*c;
    if (a <= 0.0 || disc < 0.0) return false;
    const double roots[2] = {(-b-std::sqrt(disc))/(2.0*a), (-b+std::sqrt(disc))/(2.0*a)};
    double best = forward ? std::numeric_limits<double>::infinity()
                          : -std::numeric_limits<double>::infinity();
    auto consider = [&](double t) {
        const double candidateZr = zr + t*dzr;
        if (candidateZr < -1.e-6 || candidateZr > kMoundHeightMm + 1.e-6) return;
        if (forward ? (t > 1.e-6 && t < best) : (t < -1.e-6 && t > best)) best = t;
    };
    for (double t : roots) {
        consider(t);
    }
    if (std::abs(dzr) > 1.e-6) {
        const double tBase = -zr / dzr;
        const double x = p.X() + tBase*d.X(), y = p.Y() + tBase*d.Y();
        if (x*x + y*y <= r2 + 1.e-6) consider(tBase);
    }
    if (!std::isfinite(best)) return false;
    intersection = p + best*d;
    return true;
}
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) std::cerr << "Usage: mpirun -np <cores> " << argv[0] << " <input_root_file> [skim_factor] [--ground-incline-deg D]" << std::endl;
        MPI_Finalize();
        return 1;
    }

    std::string inputFile = argv[1];

    // --- NEW: Parse the optional skim factor ---
    int skimFactor = 1; // Default is 1 (process everything)
    if (argc >= 3 && argv[2][0] != '-') {
        skimFactor = std::stoi(argv[2]);
        if (skimFactor < 1) skimFactor = 1; // Prevent division by zero
    }
    double groundInclineDeg = 0.0;
    for (int arg = 2; arg < argc; ++arg) {
        if (std::string(argv[arg]) == "--ground-incline-deg") {
            if (++arg >= argc) { if (rank == 0) std::cerr << "--ground-incline-deg requires degrees." << std::endl; MPI_Finalize(); return 1; }
            groundInclineDeg = std::stod(argv[arg]);
        }
    }
    if (!std::isfinite(groundInclineDeg) || std::abs(groundInclineDeg) >= 30.0) { if (rank == 0) std::cerr << "Invalid ground inclination." << std::endl; MPI_Finalize(); return 1; }
    const double groundSlope = std::tan(groundInclineDeg * 3.14159265358979323846 / 180.0);

    if (rank == 0) {
        std::cout << "\n--- MPI Standalone EMML Tomography Reconstruction ---\n" << std::endl;
        std::cout << "Reading from: " << inputFile << std::endl;
        if (skimFactor > 1) {
            std::cout << "SKIMMING ENABLED: Processing 1 out of every " << skimFactor << " events." << std::endl;
        }
    }


    MoundTomographyEMML_MPI emml("EMML_TomographyResults");

    emml.SetGeometryParameters(54, 54, 29, 1000.0, 1000.0, 1000.0, -8000.0);
    emml.SetMoundParameters(kMoundRadiusMm, kMoundHeightMm, groundSlope);
    emml.SetAlgorithmParameters(20, 1.0e-9, 100.0);

    TFile* fin = TFile::Open(inputFile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        if (rank == 0) std::cerr << "Error: Could not open " << inputFile << std::endl;
        MPI_Finalize();
        return 1;
    }

    TTree* tree = (TTree*)fin->Get("MuonHits");
    if (!tree) {
        std::cerr << "Rank " << rank
                  << ": Error: Could not find 'MuonHits' tree in ROOT file."
                  << std::endl;
        fin->Close();
        MPI_Finalize();
        return 1;
    }

    int totalEventsInTree = tree->GetEntries();
    
    int eventsPerCore = totalEventsInTree / size;
    int startEvent = rank * eventsPerCore;
    int endEvent = (rank == size - 1) ? totalEventsInTree : startEvent + eventsPerCore;
    
    if (rank == 0){
	    std::cout << "Total events in ROOT file: " << totalEventsInTree << std::endl;
        std::cout << "Each core will load and trace approx " << eventsPerCore << " events." << std::endl;
    }
    double in1_x, in1_y, in1_z;
    double in2_x, in2_y, in2_z;
    double in2_px, in2_py, in2_pz;
    double out1_x, out1_y, out1_z;
    double out2_x, out2_y, out2_z;
    
    tree->SetBranchAddress("In1_X", &in1_x);
    tree->SetBranchAddress("In1_Y", &in1_y);
    tree->SetBranchAddress("In1_Z", &in1_z);
    
    tree->SetBranchAddress("In2_X", &in2_x);
    tree->SetBranchAddress("In2_Y", &in2_y);
    tree->SetBranchAddress("In2_Z", &in2_z);
    
    tree->SetBranchAddress("Out1_X", &out1_x);
    tree->SetBranchAddress("Out1_Y", &out1_y);
    tree->SetBranchAddress("Out1_Z", &out1_z);
    
    tree->SetBranchAddress("Out2_X", &out2_x);
    tree->SetBranchAddress("Out2_Y", &out2_y);
    tree->SetBranchAddress("Out2_Z", &out2_z);

    tree->SetBranchAddress("In2_Px", &in2_px);
    tree->SetBranchAddress("In2_Py", &in2_py);
    tree->SetBranchAddress("In2_Pz", &in2_pz);

    std::cout << "Loading " << tree->GetEntries() << " events..." << std::endl;

    for (int i = startEvent; i < endEvent; ++i) {
        if (i % skimFactor != 0) continue;

        tree->GetEntry(i);

        TVector3 pIn1(in1_x, in1_y, in1_z);
        TVector3 pIn2(in2_x, in2_y, in2_z);
        TVector3 pOut1(out1_x, out1_y, out1_z);
        TVector3 pOut2(out2_x, out2_y, out2_z);

        MuonScatteringEvent ev;
        ev.EventId = i;
        
        ev.EntryPoint = pIn2; 
        ev.ExitPoint = pOut1; 
        
        ev.EntryDir = (pIn2 - pIn1).Unit(); 
        ev.ExitDir = (pOut2 - pOut1).Unit();  

	    
        if (!IntersectSlopedMound(pIn2, ev.EntryDir, groundSlope, true, ev.EntryPoint)
            || !IntersectSlopedMound(pOut1, ev.ExitDir, groundSlope, false, ev.ExitPoint)) {
            continue;
        }
 
        double p_magnitude = std::sqrt((in2_px * in2_px) + (in2_py * in2_py) + (in2_pz * in2_pz));
	ev.p = p_magnitude;
        
        ev.ScatteringAngle = ev.EntryDir.Angle(ev.ExitDir); 
        
        emml.AddEvent(ev);
    }
    
    fin->Close();

    // Run processing
    emml.TracePaths(); 
    emml.RunEMML();
    
    if (rank == 0) {
        // 1. Save EVERYTHING to the ROOT file first
	// min Log, max Log, min number of events
        emml.SaveResults(-20.0, 0.0, 1); 
    }
    MPI_Finalize();
    return 0;
}
