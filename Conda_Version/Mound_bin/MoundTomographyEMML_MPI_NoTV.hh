#ifndef MoundTomographyEMML_MPI_NoTV_h
#define MoundTomographyEMML_MPI_NoTV_h 1

#include <vector>
#include <string>
#include <map>
#include "TVector3.h"
#include "TH3D.h"
#include "TFile.h"

struct VoxelWeight {
    double Wx = 0.0;
    double Wy = 0.0;
    double Wxy = 0.0;
};

struct MuonScatteringEvent {
    int EventId;
    double pr; 
    double p; 
    TVector3 EntryPoint;
    TVector3 EntryDir;
    TVector3 ExitPoint;
    TVector3 ExitDir;
    double ScatteringAngle;   
    
    double DeltaThetaX;
    double DeltaThetaY;
    double DeltaX;
    double DeltaY;
    
    std::map<int, VoxelWeight> VoxelWeights; 
};

class MoundTomographyEMML_MPI {
public:
    MoundTomographyEMML_MPI(const std::string& outputFilename);
    ~MoundTomographyEMML_MPI();

    void SetGeometryParameters(int nx, int ny, int nz, double vx, double vy, double vz,
                               double zMin = 0.0);
    void SetMoundParameters(double radius, double height,
                            double groundSlope = 0.0, double baseZ = 0.0);
    void SetAlgorithmParameters(int iterations, double initLambda, double stepSize);
    void SetAngleOnly(bool enabled);
    // A value of zero retains every post-acceptance event.
    void SetMaxAcceptedMuons(long long maxMuonHits);
    
    void AddEvent(const MuonScatteringEvent& event);
    void TracePaths(); 
    void GenerateCentralCavityClosure(double soilLambda, double cavityScale,
                                      unsigned int seed);
    void PrintResidualDiagnostics(double referenceLambda) const;
    void RunEMML();
    void SaveResults(double signalMinThreshold = -10.0, double signalMaxThreshold = -5.0, int eventThreshold = 5);

private:
    std::string fOutputFilename;
    TFile* fOutputFile;
    TH3D* fRecMap3D;
    TH3D* fLogRecMap3D;
    TH3D* fScaledRecMap3D;
    TH3D* fEventMap3D;
    TH1D* fDeltaTheta;
    double fEx;
    double fEy;
    double fExy;
    bool fAngleOnly;

    int fNVoxelX, fNVoxelY, fNVoxelZ;
    double fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ;
    double fGridZMin;
    int fTotalVoxels;
    double fMoundRadius;
    double fMoundHeight;
    double fMoundBaseZ;
    double fGroundSlope;

    int fIterations;
    double fInitialLambda;
    double fstepSize;
    long long fMaxAcceptedMuons;
    std::vector<double> fLambda; 
    std::vector<int> fEventsPerVoxel;
    std::vector<MuonScatteringEvent> fEvents;

    int GetVoxelId(double x, double y, double z);
    bool IsInsideMound(double x, double y, double z) const;
};

#endif
