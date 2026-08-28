#include "MoundTomographyEMML_MPI_TV.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <mpi.h>
#include "TFile.h"
#include "TH3D.h"
#include "TVector3.h"

MoundTomographyEMML_MPI::MoundTomographyEMML_MPI(const std::string& outputFilename) 
    : fOutputFilename(outputFilename), fOutputFile(nullptr), fRecMap3D(nullptr), fLogRecMap3D(nullptr), fEventMap3D(nullptr), fDeltaTheta(nullptr),
      fIterations(10), fInitialLambda(6.11e-9), fstepSize(25.0),
      fMoundRadius(26500.0), fMoundHeight(12700.0), fGridZMin(0.0),
      fGroundSlope(0.0) {
    
    // Detector baseline resolutions. 
    // These are mathematically required to prevent the inverted covariance matrix (C^-1)
    // from exploding to infinity for short or low-density tracks.
    fEx = 1.0e-6; // rad^2
    fEy = 1.0;    // cm^2
    fExy = 0.0;   // rad*cm
}

void MoundTomographyEMML_MPI::SetMoundParameters(double radius, double height, double groundSlope) {
    fMoundRadius = radius;
    fMoundHeight = height;
    fGroundSlope = groundSlope;
}

MoundTomographyEMML_MPI::~MoundTomographyEMML_MPI() {
    if (fOutputFile) {
        fOutputFile->Close();
        delete fOutputFile;
    }
}

void MoundTomographyEMML_MPI::SetGeometryParameters(int nx, int ny, int nz, double vx, double vy, double vz, double zMin) {
    fNVoxelX = nx; fNVoxelY = ny; fNVoxelZ = nz;
    fVoxelSizeX = vx; fVoxelSizeY = vy; fVoxelSizeZ = vz;
    fGridZMin = zMin;
    fTotalVoxels = nx * ny * nz;

    fLambda.assign(fTotalVoxels, fInitialLambda);
    fEventsPerVoxel.assign(fTotalVoxels, 0);
    
    // Pre-initialize air voxels to the noise floor to prevent TV boundary ringing
    for (int iz = 0; iz < fNVoxelZ; ++iz) {
        for (int iy = 0; iy < fNVoxelY; ++iy) {
            for (int ix = 0; ix < fNVoxelX; ++ix) {
                // Calculate voxel center coordinates
                double x = (ix + 0.5) * fVoxelSizeX - (fNVoxelX * fVoxelSizeX / 2.0);
                double y = (iy + 0.5) * fVoxelSizeY - (fNVoxelY * fVoxelSizeY / 2.0);
                double z = fGridZMin + (iz + 0.5) * fVoxelSizeZ;
                
                if (!IsInsideMound(x, y, z)) {
                    int idx = ix + iy * fNVoxelX + iz * fNVoxelX * fNVoxelY;
                    fLambda[idx] = 1e-12; // Air density floor
                }
            }
        }
    }
}

void MoundTomographyEMML_MPI::SetAlgorithmParameters(int iterations, double initialLambda, double stepSize) {
    fIterations = iterations;
    fInitialLambda = initialLambda;
    fstepSize = stepSize;
    
    // Re-initialize the density map in case the initial guess changed
    if (fTotalVoxels > 0) {
        // Safely re-apply the geometry and air mask
        SetGeometryParameters(fNVoxelX, fNVoxelY, fNVoxelZ, fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ, fGridZMin);
    }
}

int MoundTomographyEMML_MPI::GetVoxelId(double x, double y, double z) {
    // Centers the coordinate grid on X and Y; Z is in global Geant4 coordinates.
    int ix = std::floor((x - (-fNVoxelX * fVoxelSizeX / 2.0)) / fVoxelSizeX);
    int iy = std::floor((y - (-fNVoxelY * fVoxelSizeY / 2.0)) / fVoxelSizeY);
    int iz = std::floor((z - fGridZMin) / fVoxelSizeZ);

    if (ix >= 0 && ix < fNVoxelX && iy >= 0 && iy < fNVoxelY && iz >= 0 && iz < fNVoxelZ) {
        return ix + iy * fNVoxelX + iz * fNVoxelX * fNVoxelY;
    }
    return -1;
}

bool MoundTomographyEMML_MPI::IsInsideMound(double x, double y, double z) const {
    const double zRelativeToGround = z + fGroundSlope * y;
    if (zRelativeToGround < 0.0 || zRelativeToGround > fMoundHeight) return false;

    const double radialTerm = (x * x + y * y) / (fMoundRadius * fMoundRadius);
    const double verticalTerm = (zRelativeToGround * zRelativeToGround) /
                                (fMoundHeight * fMoundHeight);
    return radialTerm + verticalTerm <= 1.0;
}

void MoundTomographyEMML_MPI::AddEvent(const MuonScatteringEvent& ev) {
    MuonScatteringEvent event = ev;

    // 1. Define Local 3D Track Coordinate System 
    TVector3 u = event.EntryDir.Unit();
    TVector3 z_axis(0.0, 0.0, 1.0);
    TVector3 v = u.Cross(z_axis);

    if (v.Mag() < 1e-6) {
        v = u.Cross(TVector3(1.0, 0.0, 0.0));
    }
    v = v.Unit();
    TVector3 w = u.Cross(v).Unit();

    // 2. Local Angular Deviation
    event.DeltaThetaX = std::asin(event.ExitDir.Dot(v));
    event.DeltaThetaY = std::asin(event.ExitDir.Dot(w));

    // 3. Calculate Local Spatial Displacement 
    TVector3 deltaP = event.ExitPoint - event.EntryPoint;
    double L = deltaP.Dot(u); 
    
    TVector3 expectedPoint = event.EntryPoint + (L * u);
    TVector3 spatialDisplacement = event.ExitPoint - expectedPoint;

    // Split the deviation onto our orthogonal planes and convert to cm
    event.DeltaX = spatialDisplacement.Dot(v) / 10.0;
    event.DeltaY = spatialDisplacement.Dot(w) / 10.0;

    fEvents.push_back(event);
}

void MoundTomographyEMML_MPI::TracePaths() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Limit the number of accepted MuonScatteringEvents, not their EventId.
    // The exclusive prefix sum gives each rank its position in the global
    // event stream, so the cap remains exactly 2,500,000 with any MPI size.
    constexpr long long kMaxMuonHits = 2500000;
    const long long localEventCount = static_cast<long long>(fEvents.size());
    long long precedingEventCount = 0;
    MPI_Exscan(&localEventCount, &precedingEventCount, 1, MPI_LONG_LONG,
               MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) precedingEventCount = 0;

    if (precedingEventCount >= kMaxMuonHits) {
        fEvents.clear();
    } else if (precedingEventCount + localEventCount > kMaxMuonHits) {
        fEvents.resize(static_cast<std::size_t>(kMaxMuonHits - precedingEventCount));
    }

    const long long retainedLocalCount = static_cast<long long>(fEvents.size());
    long long retainedGlobalCount = 0;
    MPI_Allreduce(&retainedLocalCount, &retainedGlobalCount, 1, MPI_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Tracing " << retainedGlobalCount
                  << " events mathematically..." << std::endl;
    }

    const double xMin = -fNVoxelX * fVoxelSizeX / 2.0;
    const double yMin = -fNVoxelY * fVoxelSizeY / 2.0;
    const double zMin = fGridZMin;
    const double xMax = -xMin;
    const double yMax = -yMin;
    const double zMax = fGridZMin + fNVoxelZ * fVoxelSizeZ;
    const double eps = 1.0e-9;

    for (auto& event : fEvents) {
        event.VoxelWeights.clear();
        TVector3 direction = (event.ExitPoint - event.EntryPoint).Unit();
        double totalLength = (event.ExitPoint - event.EntryPoint).Mag();
        if (totalLength <= eps) continue;

        double tEnter = 0.0;
        double tExit = totalLength;
        const double p[3] = {event.EntryPoint.X(), event.EntryPoint.Y(), event.EntryPoint.Z()};
        const double d[3] = {direction.X(), direction.Y(), direction.Z()};
        const double lower[3] = {xMin, yMin, zMin};
        const double upper[3] = {xMax, yMax, zMax};
        bool intersects = true;
        for (int axis = 0; axis < 3; ++axis) {
            if (std::abs(d[axis]) < eps) {
                if (p[axis] < lower[axis] || p[axis] > upper[axis]) intersects = false;
                continue;
            }
            double a = (lower[axis] - p[axis]) / d[axis];
            double b = (upper[axis] - p[axis]) / d[axis];
            if (a > b) std::swap(a, b);
            tEnter = std::max(tEnter, a);
            tExit = std::min(tExit, b);
        }
        if (!intersects || tExit <= tEnter + eps) continue;

        TVector3 currentPos = event.EntryPoint + direction * (tEnter + eps);
        int voxelId = GetVoxelId(currentPos.X(), currentPos.Y(), currentPos.Z());
        if (voxelId < 0) continue;

        int ix = voxelId % fNVoxelX;
        int iy = (voxelId / fNVoxelX) % fNVoxelY;
        int iz = voxelId / (fNVoxelX * fNVoxelY);
        int index[3] = {ix, iy, iz};
        const int nVoxel[3] = {fNVoxelX, fNVoxelY, fNVoxelZ};
        const double voxelSize[3] = {fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ};

        int step[3];
        double tNext[3];
        double tDelta[3];
        for (int axis = 0; axis < 3; ++axis) {
            if (std::abs(d[axis]) < eps) {
                step[axis] = 0;
                tNext[axis] = std::numeric_limits<double>::infinity();
                tDelta[axis] = std::numeric_limits<double>::infinity();
            } else {
                step[axis] = d[axis] > 0.0 ? 1 : -1;
                double nextBoundary = lower[axis] + (index[axis] + (step[axis] > 0 ? 1 : 0)) * voxelSize[axis];
                tNext[axis] = (nextBoundary - p[axis]) / d[axis];
                tDelta[axis] = voxelSize[axis] / std::abs(d[axis]);
            }
        }

        double t = tEnter;
        double accumulated_length_cm = 0.0; // Track the total length inside the grid
        
        while (t < tExit - eps) {
            const double nextCrossing = std::min(tNext[0], std::min(tNext[1], tNext[2]));
            const double tStop = std::min(nextCrossing, tExit);
            const double segmentMm = tStop - t;
            if (segmentMm > eps) {
                const int id = index[0] + index[1] * fNVoxelX + index[2] * fNVoxelX * fNVoxelY;
                const double lengthCm = segmentMm / 10.0;
                
                accumulated_length_cm += lengthCm; // Accumulate the chord length
                
                const double remainingCm = (totalLength - tStop) / 10.0;
                event.VoxelWeights[id].Wx += lengthCm;
                event.VoxelWeights[id].Wxy += lengthCm * remainingCm + 0.5 * lengthCm * lengthCm;
                event.VoxelWeights[id].Wy += lengthCm * remainingCm * remainingCm
                                            + lengthCm * lengthCm * remainingCm
                                            + lengthCm * lengthCm * lengthCm / 3.0;
                if (event.VoxelWeights[id].Wx == lengthCm) ++fEventsPerVoxel[id];
            }
            t = tStop;
            for (int axis = 0; axis < 3; ++axis) {
                if (std::abs(tNext[axis] - nextCrossing) < eps) {
                    index[axis] += step[axis];
                    tNext[axis] += tDelta[axis];
                }
            }
            if (index[0] < 0 || index[0] >= nVoxel[0] || index[1] < 0 || index[1] >= nVoxel[1] || index[2] < 0 || index[2] >= nVoxel[2]) break;
        }
        
        // MINIMUM CHORD CUT: Drop tracks that merely graze the boundaries
        // 100.0 cm (1 meter) threshold filters out skimming tracks causing boundary spikes.
        if (accumulated_length_cm < 100.0) {
            // Undo the event counts for the voxels touched by this grazing track
            for (const auto& weight_pair : event.VoxelWeights) {
                --fEventsPerVoxel[weight_pair.first]; 
            }
            // Clear the weights so RunEMML naturally skips this event
            event.VoxelWeights.clear();
        }
    }
}

void MoundTomographyEMML_MPI::RunEMML() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> local_numerator(fTotalVoxels, 0.0);
    std::vector<double> local_denominator(fTotalVoxels, 0.0);
    std::vector<double> global_numerator(fTotalVoxels, 0.0);
    std::vector<double> global_denominator(fTotalVoxels, 0.0);

    std::vector<int> global_events_per_voxel(fTotalVoxels, 0);
    MPI_Allreduce(fEventsPerVoxel.data(), global_events_per_voxel.data(), fTotalVoxels,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    fEventsPerVoxel.swap(global_events_per_voxel);

    if (rank == 0) std::cout << "Starting MAP-EM iterations with TV Regularization..." << std::endl;

    for (int iter = 0; iter < fIterations; ++iter) {
        std::fill(local_numerator.begin(), local_numerator.end(), 0.0);
        std::fill(local_denominator.begin(), local_denominator.end(), 0.0);

        for (const auto& event : fEvents) {
            if (event.VoxelWeights.empty()) continue; // Skip dropped grazing tracks

            constexpr double kReferenceMomentumMeV = 10000.0;
            const double momentumScale = event.p > 0.0
                ? std::pow(kReferenceMomentumMeV / event.p, 2)
                : 1.0;
                
            // 1. Build the 2x2 Expected Covariance Matrix
            double C_angle = 0.0; 
            double C_pos   = 0.0; 
            double C_cov   = 0.0; 

            for (const auto& weight_pair : event.VoxelWeights) {
                int voxelId = weight_pair.first;
                const VoxelWeight& weights = weight_pair.second;

                C_angle += momentumScale * weights.Wx  * fLambda[voxelId];
                C_pos   += momentumScale * weights.Wy  * fLambda[voxelId];
                C_cov   += momentumScale * weights.Wxy * fLambda[voxelId];
            }

            // 2. Add Detector Baseline Noise
            C_angle += fEx;  
            C_pos   += fEy;  
            C_cov   += fExy; 

            // 3. Invert the 2x2 Matrix
            double det = (C_angle * C_pos) - (C_cov * C_cov);
            if (!std::isfinite(det) || det <= 0.0) continue;
            
            double inv_C_angle =  C_pos / det;
            double inv_C_pos   =  C_angle / det;
            double inv_C_cov   = -C_cov / det;

            for (const auto& weight_pair : event.VoxelWeights) {
                int voxelId = weight_pair.first;
                const VoxelWeight& weights = weight_pair.second;

                const double wx = momentumScale * weights.Wx;
                const double wy = momentumScale * weights.Wy;
                const double wxy = momentumScale * weights.Wxy;
                const double q00 = inv_C_angle * inv_C_angle * wx
                                  + 2.0 * inv_C_angle * inv_C_cov * wxy
                                  + inv_C_cov * inv_C_cov * wy;
                const double q01 = inv_C_angle * inv_C_cov * wx
                                  + (inv_C_angle * inv_C_pos + inv_C_cov * inv_C_cov) * wxy
                                  + inv_C_cov * inv_C_pos * wy;
                const double q11 = inv_C_cov * inv_C_cov * wx
                                  + 2.0 * inv_C_cov * inv_C_pos * wxy
                                  + inv_C_pos * inv_C_pos * wy;
                const double trace = inv_C_angle * wx
                                   + 2.0 * inv_C_cov * wxy
                                   + inv_C_pos * wy;

                const double scoreX = event.DeltaThetaX * event.DeltaThetaX * q00
                                    + 2.0 * event.DeltaThetaX * event.DeltaX * q01
                                    + event.DeltaX * event.DeltaX * q11;
                const double scoreY = event.DeltaThetaY * event.DeltaThetaY * q00
                                    + 2.0 * event.DeltaThetaY * event.DeltaY * q01
                                    + event.DeltaY * event.DeltaY * q11;

                local_numerator[voxelId] += scoreX + scoreY;
                local_denominator[voxelId] += 2.0 * trace;
            }
        }

        // MPI Synchronization across cores
        MPI_Allreduce(local_numerator.data(), global_numerator.data(), fTotalVoxels, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local_denominator.data(), global_denominator.data(), fTotalVoxels, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // -------------------------------------------------------------
        // MAP-EM Maximization Update Step with 3D TV Regularization
        // -------------------------------------------------------------
        // beta acts as a relative strength parameter (e.g. 0.05 = 5% smoothing pressure)
        const double beta = 0.05;   
        // delta is scaled closer to true physical lambda values (~10^-9)
        const double delta = 1e-11; 
        const double deltaSq = delta * delta;
        std::vector<double> lambda_old = fLambda;

        for (int iz = 0; iz < fNVoxelZ; ++iz) {
            for (int iy = 0; iy < fNVoxelY; ++iy) {
                for (int ix = 0; ix < fNVoxelX; ++ix) {
                    
                    int i = ix + iy * fNVoxelX + iz * fNVoxelX * fNVoxelY;
                    if (global_denominator[i] <= 0) continue;
                    // CRITICAL FIX: Do not update air voxels!
                    // Let them remain at their 1e-12 floor without participating in EMML.
                    double cx = (ix + 0.5) * fVoxelSizeX - (fNVoxelX * fVoxelSizeX / 2.0);
                    double cy = (iy + 0.5) * fVoxelSizeY - (fNVoxelY * fVoxelSizeY / 2.0);
                    double cz = fGridZMin + (iz + 0.5) * fVoxelSizeZ;
                    if (!IsInsideMound(cx, cy, cz)) continue;
                    double val_j = lambda_old[i];
                    double tv_gradient = 0.0;

                    // 6-Connectivity Neighbor Offsets in 3D
                    int neighbor_indices[6][3] = {
                        {ix - 1, iy, iz}, {ix + 1, iy, iz}, // Left, Right
                        {ix, iy - 1, iz}, {ix, iy + 1, iz}, // Front, Back
                        {ix, iy, iz - 1}, {ix, iy, iz + 1}  // Down, Up
                    };

                    for (int n = 0; n < 6; ++n) {
                        int nx = neighbor_indices[n][0];
                        int ny = neighbor_indices[n][1];
                        int nz = neighbor_indices[n][2];

                        // Boundary condition check
                        if (nx >= 0 && nx < fNVoxelX && ny >= 0 && ny < fNVoxelY && nz >= 0 && nz < fNVoxelZ) {
                            int neighbor_idx = nx + ny * fNVoxelX + nz * fNVoxelX * fNVoxelY;

                            // Disconnect TV smoothing across the physical boundary!
                            // Only smooth dirt against dirt. If a neighbor is outside the mound (air),
                            // ignore it so it doesn't drag the dirt density down to 1e-12.
                            double nx_coord = (nx + 0.5) * fVoxelSizeX - (fNVoxelX * fVoxelSizeX / 2.0);
                            double ny_coord = (ny + 0.5) * fVoxelSizeY - (fNVoxelY * fVoxelSizeY / 2.0);
                            double nz_coord = (nz + 0.5) * fVoxelSizeZ;
                            
                            if (!IsInsideMound(nx_coord, ny_coord, nz_coord)) continue;
			    if (global_denominator[neighbor_idx] <= 0) continue;


                            double val_k = lambda_old[neighbor_idx];
                            double diff = val_j - val_k;
                            
                            // Accumulate the smoothed derivative component
                            tv_gradient += diff / std::sqrt(diff * diff + deltaSq);
                        }
                    }

                    // Preconditioned Penalty:
                    double penalty = beta * global_denominator[i] * (tv_gradient / 6.0);

                    // Calculate regularized denominator
                    double regularized_denominator = global_denominator[i] + penalty;

                    // OSL stability fallback: If the penalty is aggressively negative, 
                    // prevent the denominator from dropping below 10% of its physical data value.
                    if (regularized_denominator < 0.1 * global_denominator[i]) {
                        regularized_denominator = 0.1 * global_denominator[i];
                    }

                    double update_factor = global_numerator[i] / regularized_denominator;
                    if (std::isfinite(update_factor) && update_factor >= 0.0) {
                        fLambda[i] *= update_factor;
                    }

                    // Enforce a strict positive physical density floor.
//                    if (fLambda[i] < 1e-12) {
//                        fLambda[i] = 1e-12;
//                    }
                }
            }
        }
        
        if (rank == 0) {
            std::cout << "Completed Iteration " << iter + 1 << "/" << fIterations << std::endl;
        }
    }
}

void MoundTomographyEMML_MPI::SaveResults(double minThreshold, double maxThreshold, int minEvents) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) return;

    double xMin = -fNVoxelX * fVoxelSizeX / 2.0;
    double xMax =  fNVoxelX * fVoxelSizeX / 2.0;
    double yMin = -fNVoxelY * fVoxelSizeY / 2.0;
    double yMax =  fNVoxelY * fVoxelSizeY / 2.0;
    double zMin = fGridZMin;
    double zMax = fGridZMin + fNVoxelZ * fVoxelSizeZ;

    fRecMap3D = new TH3D("RecMap3D", "Reconstructed Density Map", fNVoxelX, xMin, xMax, fNVoxelY, yMin, yMax, fNVoxelZ, zMin, zMax);
    fLogRecMap3D = new TH3D("LogRecMap3D", "Log(Lambda) Map", fNVoxelX, xMin, xMax, fNVoxelY, yMin, yMax, fNVoxelZ, zMin, zMax);
    fEventMap3D = new TH3D("EventMap3D", "Event Count Map", fNVoxelX, xMin, xMax, fNVoxelY, yMin, yMax, fNVoxelZ, zMin, zMax);

    FILE* fout_vtk = fopen((fOutputFilename + ".vtk").c_str(), "w");
    if (fout_vtk) {
        fprintf(fout_vtk, "# vtk DataFile Version 2.0\n");
        fprintf(fout_vtk, "Tomography Rectilinear Grid\n");    
        fprintf(fout_vtk, "ASCII\n");    
        fprintf(fout_vtk, "DATASET RECTILINEAR_GRID\n");    
        fprintf(fout_vtk, "DIMENSIONS %d %d %d\n", fNVoxelX, fNVoxelY, fNVoxelZ);    
        
        fprintf(fout_vtk, "X_COORDINATES %d float\n", fNVoxelX);    
        for(int i=0; i<fNVoxelX; i++) fprintf(fout_vtk, "%f\n", fRecMap3D->GetXaxis()->GetBinLowEdge(i+1));        
        
        fprintf(fout_vtk, "Y_COORDINATES %d float\n", fNVoxelY);    
        for(int j=0; j<fNVoxelY; j++) fprintf(fout_vtk, "%f\n", fRecMap3D->GetYaxis()->GetBinLowEdge(j+1));            
        
        fprintf(fout_vtk, "Z_COORDINATES %d float\n", fNVoxelZ);
        for(int k=0; k<fNVoxelZ; k++) fprintf(fout_vtk, "%f\n", fRecMap3D->GetZaxis()->GetBinLowEdge(k+1));
        
        fprintf(fout_vtk, "POINT_DATA %d\n", fTotalVoxels);
        fprintf(fout_vtk, "SCALARS Density float\n");
        fprintf(fout_vtk, "LOOKUP_TABLE default\n");
        
        for(int k=1; k<=fNVoxelZ; k++) {
            for(int j=1; j<=fNVoxelY; j++) {
                for(int i=1; i<=fNVoxelX; i++) {
                    int idx = (i-1) + (j-1)*fNVoxelX + (k-1)*fNVoxelX*fNVoxelY;
                    constexpr double kMaskedLogDensity = -15;
                    constexpr double kMaskedDensity = 1e-15;
                    double lambda = fLambda[idx];
                    int events = fEventsPerVoxel[idx];
                    const double x = fRecMap3D->GetXaxis()->GetBinCenter(i);
                    const double y = fRecMap3D->GetYaxis()->GetBinCenter(j);
                    const double z = fRecMap3D->GetZaxis()->GetBinCenter(k);
                    const bool outsideMound = !IsInsideMound(x, y, z);
                    const double logLambda = outsideMound
                        ? kMaskedLogDensity
                        : std::log10(lambda);
                    const double Lambda = outsideMound
                        ? kMaskedDensity
                        : lambda;


                    fRecMap3D->SetBinContent(i, j, k, lambda);
                    // Keep the ROOT and VTK log fields consistent: voxels
                    // outside the physical mound are a mask, not reconstructed air.
                    fLogRecMap3D->SetBinContent(i, j, k,
                                                 outsideMound ? kMaskedLogDensity : logLambda);
                    fEventMap3D->SetBinContent(i, j, k, events);
                    // Clamp limits for visual exports
                    if (outsideMound) {
			//fprintf(fout_vtk, "%g\n", kMaskedLogDensity);
			fprintf(fout_vtk, "%g\n", -1.0000);
                    }
		    else if (events < minEvents) {
			//fprintf(fout_vtk, "%g\n", kMaskedLogDensity);
			fprintf(fout_vtk, "%g\n", 0.0000);
                    }
		    else if (logLambda < minThreshold || logLambda > maxThreshold) {
			//fprintf(fout_vtk, "%g\n", kMaskedLogDensity);
			fprintf(fout_vtk, "%g\n", kMaskedDensity);
                    }
		    else {
                        //fprintf(fout_vtk, "%e\n", logLambda);
                        fprintf(fout_vtk, "%e\n", Lambda);
		    }
                }
            }
        }
        fclose(fout_vtk);
    }
    
    fOutputFile = new TFile((fOutputFilename + ".root").c_str(), "RECREATE");
    fRecMap3D->Write();
    fLogRecMap3D->Write();
    fEventMap3D->Write();
    fOutputFile->Close();
}
