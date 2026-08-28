#include "MoundTomographyEMML_MPI_NoTV.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH3D.h"
#include "TH1D.h"
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <mpi.h>

namespace {
constexpr double kMoundRadiusMm = 26500.0;
constexpr double kHorizontalMoundHeightMm = 12700.0;
constexpr double kInclinedMoundClearanceMm = 15000.0;
constexpr double kMoundGroundSeparationMm = 0.001;

double MoundBaseZMm(double slope) {
    return std::abs(slope) < 1.e-15 ? 0.0
                                    : -std::abs(slope) * kMoundRadiusMm;
}

double MoundHeightMm(double slope) {
    return std::abs(slope) < 1.e-15
        ? kHorizontalMoundHeightMm
        : kInclinedMoundClearanceMm
          + 2.0 * std::abs(slope) * kMoundRadiusMm;
}

double MoundGroundSeparationMm(double slope) {
    return std::abs(slope) < 1.e-15 ? 0.0 : kMoundGroundSeparationMm;
}

// Intersect the same vertical analytic ellipsoid clipped above
// z=-slope*y that is constructed in DetectorConstruction.cc. `forward`
// chooses the nearest positive (entry) or negative (exit) ray parameter.
bool IntersectMound(const TVector3& p, const TVector3& d, double slope,
                    bool forward, TVector3& intersection) {
    const double moundBaseZ = MoundBaseZMm(slope);
    const double moundHeight = MoundHeightMm(slope);
    const double zRel = p.Z() - moundBaseZ;
    const double dzRel = d.Z();
    const double r2 = kMoundRadiusMm * kMoundRadiusMm;
    const double h2 = moundHeight * moundHeight;
    const double a = (d.X()*d.X() + d.Y()*d.Y()) / r2 + dzRel*dzRel / h2;
    const double b = 2.0 * ((p.X()*d.X() + p.Y()*d.Y()) / r2 + zRel*dzRel / h2);
    const double c = (p.X()*p.X() + p.Y()*p.Y()) / r2 + zRel*zRel / h2 - 1.0;
    const double disc = b*b - 4.0*a*c;
    if (a <= 0.0 || disc < 0.0) return false;

    const double roots[2] = {(-b - std::sqrt(disc)) / (2.0*a),
                             (-b + std::sqrt(disc)) / (2.0*a)};
    const double eps = 1.e-6;
    double chosen = forward ? std::numeric_limits<double>::infinity()
                            : -std::numeric_limits<double>::infinity();
    auto consider = [&](double t) {
        const double candidateZRel = zRel + t * dzRel;
        const double candidateGroundRel =
            p.Z() + t*d.Z() + slope*(p.Y() + t*d.Y())
            - MoundGroundSeparationMm(slope);
        if (candidateZRel < -eps || candidateZRel > moundHeight + eps
            || candidateGroundRel < -eps) {
            return;
        }
        if (forward ? (t > eps && t < chosen) : (t < -eps && t > chosen)) chosen = t;
    };
    for (double t : roots) {
        consider(t);
    }
    // The Boolean mound is closed by the inclined ground plane.
    const double groundRel =
        p.Z() + slope*p.Y() - MoundGroundSeparationMm(slope);
    const double dGroundRel = d.Z() + slope*d.Y();
    if (std::abs(dGroundRel) > eps) {
        const double tBase = -groundRel / dGroundRel;
        const double x = p.X() + tBase*d.X();
        const double y = p.Y() + tBase*d.Y();
        const double candidateZRel = p.Z() + tBase*d.Z() - moundBaseZ;
        const double ellipsoidValue =
            (x*x + y*y)/r2 + candidateZRel*candidateZRel/h2;
        if (ellipsoidValue <= 1.0 + eps) consider(tBase);
    }
    if (!std::isfinite(chosen)) return false;
    intersection = p + chosen * d;
    return true;
}
} // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) std::cerr << "Usage: mpirun -np <cores> " << argv[0]
                                 << " <input_root_file> [skim_factor] [--closure [seed]]"
                                 << " [--diagnostics-only] [--diagnostics-lambda value]"
                                 << " [--angle-only] [--initial-lambda value]"
                                 << " [--virtual-four-patch] [--sparse-four-pair]"
                                 << " [--sparse-eight-pair]"
                                 << " [--detector-pairs N --detector-width M --detector-height M"
                                 << " --detector-elevation M --detector-incline-deg D"
                                 << " --detector-pair-angle-deg A (repeat N times)]"
                                 << " [--min-muon-horizontal-y C]"
                                 << " [--ground-incline-deg D]"
                                 << " [--voxel-size-y M]"
                                 << " [--max-accepted-muons N]"
                                 << std::endl;
        MPI_Finalize();
        return 1;
    }

    std::string inputFile = argv[1];

    // Parse the optional skim factor and closure-test mode.
    int skimFactor = 1; // Default is 1 (process everything)
    bool closureMode = false;
    bool diagnosticsOnly = false;
    bool angleOnly = false;
    bool virtualFourPatch = false;
    bool sparseFourPair = false;
    bool sparseEightPair = false;
    bool detectorPairsSpecified = false;
    bool detectorWidthSpecified = false;
    bool detectorHeightSpecified = false;
    bool detectorElevationSpecified = false;
    bool detectorInclineSpecified = false;
    double groundInclineDeg = 0.0;
    int detectorPairs = 0;
    double detectorWidthM = 0.0;
    double detectorHeightM = 0.0;
    double detectorElevationM = 0.0;
    std::vector<double> detectorPairAnglesDeg;
    double minMuonHorizontalY = -1.0;
    double detectorInclineDeg = 0.0;
    double diagnosticsLambda = 1.0e-7;
    double initialLambda = 1.0e-9;
    double voxelSizeYM = 1.0;
    long long maxAcceptedMuons = 2500000;
    unsigned int closureSeed = 12345U;
    for (int arg = 2; arg < argc; ++arg) {
        const std::string option = argv[arg];
        if (option == "--closure") {
            closureMode = true;
            if (arg + 1 < argc && argv[arg + 1][0] != '-') {
                closureSeed = static_cast<unsigned int>(std::stoul(argv[++arg]));
            }
        } else if (option == "--diagnostics-only") {
            diagnosticsOnly = true;
        } else if (option == "--diagnostics-lambda") {
            if (arg + 1 >= argc) {
                if (rank == 0) {
                    std::cerr << "Error: --diagnostics-lambda requires a positive value."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            try {
                diagnosticsLambda = std::stod(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) {
                    std::cerr << "Error: invalid diagnostic lambda value."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            if (diagnosticsLambda <= 0.0) {
                if (rank == 0) {
                    std::cerr << "Error: diagnostic lambda must be positive."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            diagnosticsOnly = true;
        } else if (option == "--angle-only") {
            angleOnly = true;
        } else if (option == "--virtual-four-patch") {
            virtualFourPatch = true;
        } else if (option == "--sparse-four-pair") {
            sparseFourPair = true;
        } else if (option == "--sparse-eight-pair") {
            sparseEightPair = true;
        } else if (option == "--detector-pairs") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --detector-pairs requires an integer." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                detectorPairs = std::stoi(argv[++arg]);
                detectorPairsSpecified = true;
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid detector pair count." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--detector-width") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --detector-width requires metres." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                detectorWidthM = std::stod(argv[++arg]);
                detectorWidthSpecified = true;
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid detector width." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--detector-height") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --detector-height requires metres." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                detectorHeightM = std::stod(argv[++arg]);
                detectorHeightSpecified = true;
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid detector height." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--detector-elevation") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --detector-elevation requires metres." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                detectorElevationM = std::stod(argv[++arg]);
                detectorElevationSpecified = true;
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid detector elevation." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (!std::isfinite(detectorElevationM)) {
                if (rank == 0) std::cerr << "Error: detector elevation must be finite." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--detector-pair-angle-deg") {
            if (arg + 1 >= argc) {
                if (rank == 0) {
                    std::cerr << "Error: --detector-pair-angle-deg requires degrees."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            try {
                const double angleDeg = std::stod(argv[++arg]);
                if (!std::isfinite(angleDeg)) {
                    throw std::invalid_argument("non-finite detector pair angle");
                }
                detectorPairAnglesDeg.push_back(angleDeg);
            } catch (const std::exception&) {
                if (rank == 0) {
                    std::cerr << "Error: invalid detector pair angle."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--min-muon-horizontal-y") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --min-muon-horizontal-y requires a value." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                minMuonHorizontalY = std::stod(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid minimum horizontal Y direction." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (!std::isfinite(minMuonHorizontalY)
                || minMuonHorizontalY < -1.0 || minMuonHorizontalY > 1.0) {
                if (rank == 0) {
                    std::cerr << "Error: minimum horizontal Y direction must be in [-1, 1]."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--detector-incline-deg") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --detector-incline-deg requires degrees." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                detectorInclineDeg = std::stod(argv[++arg]);
                detectorInclineSpecified = true;
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid detector inclination." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (!std::isfinite(detectorInclineDeg) || std::abs(detectorInclineDeg) >= 45.0) {
                if (rank == 0) {
                    std::cerr << "Error: detector inclination must be finite and between -45 and 45 degrees."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--ground-incline-deg") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --ground-incline-deg requires degrees." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                groundInclineDeg = std::stod(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid ground inclination." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (!std::isfinite(groundInclineDeg) || std::abs(groundInclineDeg) >= 30.0) {
                if (rank == 0) std::cerr << "Error: ground inclination must be finite and between -30 and 30 degrees." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--max-accepted-muons") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --max-accepted-muons requires a non-negative integer." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                maxAcceptedMuons = std::stoll(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid maximum accepted muon count." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (maxAcceptedMuons < 0) {
                if (rank == 0) std::cerr << "Error: maximum accepted muon count cannot be negative." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--voxel-size-y") {
            if (arg + 1 >= argc) {
                if (rank == 0) std::cerr << "Error: --voxel-size-y requires metres." << std::endl;
                MPI_Finalize();
                return 1;
            }
            try {
                voxelSizeYM = std::stod(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) std::cerr << "Error: invalid Y voxel size." << std::endl;
                MPI_Finalize();
                return 1;
            }
            if (!std::isfinite(voxelSizeYM) || voxelSizeYM <= 0.0) {
                if (rank == 0) std::cerr << "Error: Y voxel size must be finite and positive." << std::endl;
                MPI_Finalize();
                return 1;
            }
        } else if (option == "--initial-lambda") {
            if (arg + 1 >= argc) {
                if (rank == 0) {
                    std::cerr << "Error: --initial-lambda requires a positive value."
                              << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            try {
                initialLambda = std::stod(argv[++arg]);
            } catch (const std::exception&) {
                if (rank == 0) {
                    std::cerr << "Error: invalid initial lambda value." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
            if (initialLambda <= 0.0) {
                if (rank == 0) {
                    std::cerr << "Error: initial lambda must be positive." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else {
            skimFactor = std::stoi(option);
            if (skimFactor < 1) skimFactor = 1;
        }
    }
    const double groundSlope = std::tan(groundInclineDeg * 3.14159265358979323846 / 180.0);

    const bool genericDetectorMode = detectorPairsSpecified || detectorWidthSpecified
                                  || detectorHeightSpecified || detectorElevationSpecified
                                  || detectorInclineSpecified
                                  || !detectorPairAnglesDeg.empty();
    if (genericDetectorMode) {
        if (!detectorPairsSpecified || !detectorWidthSpecified || !detectorHeightSpecified
            || detectorPairs < 1 || detectorWidthM <= 0.0 || detectorHeightM <= 0.0) {
            if (rank == 0) {
                std::cerr << "Error: virtual detector mode requires positive --detector-pairs, "
                          << "--detector-width, and --detector-height values." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
        if (detectorPairAnglesDeg.size()
            != static_cast<std::size_t>(detectorPairs)) {
            if (rank == 0) {
                std::cerr << "Error: --detector-pairs " << detectorPairs
                          << " requires exactly " << detectorPairs
                          << " --detector-pair-angle-deg values, but received "
                          << detectorPairAnglesDeg.size() << "." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
        // Pair axes are equivalent modulo 180 degrees because each pair has
        // two opposing stacks. Reject duplicate physical placements.
        for (std::size_t i = 0; i < detectorPairAnglesDeg.size(); ++i) {
            double canonicalI = std::fmod(detectorPairAnglesDeg[i], 180.0);
            if (canonicalI < -90.0) canonicalI += 180.0;
            if (canonicalI >= 90.0) canonicalI -= 180.0;
            for (std::size_t j = 0; j < i; ++j) {
                double canonicalJ = std::fmod(detectorPairAnglesDeg[j], 180.0);
                if (canonicalJ < -90.0) canonicalJ += 180.0;
                if (canonicalJ >= 90.0) canonicalJ -= 180.0;
                if (std::abs(canonicalI - canonicalJ) < 1.e-9) {
                    if (rank == 0) {
                        std::cerr << "Error: detector pair angles "
                                  << detectorPairAnglesDeg[j] << " and "
                                  << detectorPairAnglesDeg[i]
                                  << " degrees describe the same opposing pair."
                                  << std::endl;
                    }
                    MPI_Finalize();
                    return 1;
                }
            }
        }
        // The physical detector is a rigid cylinder normal to the ground, so
        // its virtual panel selection defaults to that same inclination.
        if (!detectorInclineSpecified) {
            detectorInclineDeg = groundInclineDeg;
        }
        if (std::abs(detectorInclineDeg - groundInclineDeg) > 1.e-9) {
            if (rank == 0) {
                std::cerr << "Error: detector inclination must equal ground "
                          << "inclination for the rigid detector cylinder."
                          << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
        // Elevation is measured along the tilted cylinder axis from its lower
        // cap in the ground plane.
        if (!detectorElevationSpecified) {
            detectorElevationM = detectorHeightM / 2.0;
        }
        if (virtualFourPatch || sparseFourPair || sparseEightPair) {
            if (rank == 0) {
                std::cerr << "Error: generic detector dimensions cannot be combined with a "
                          << "legacy virtual detector preset." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
    }

    // Match the Geant4 analytic mound envelope. Keeping the bounds on one-metre
    // edges makes runs at different inclinations overlayable.
    constexpr double kVoxelSizeXMm = 1000.0;
    constexpr double kVoxelSizeZMm = 1000.0;
    const double voxelSizeYMm = voxelSizeYM * 1000.0;
    // Use an even number of Y cells so y=0 and the room edges at y=+/-4 m
    // are voxel boundaries for sizes such as 1, 2, and 4 m. Round upward
    // until the complete 53 m mound diameter remains inside the grid.
    const double requiredVoxelCountY =
        std::ceil(2.0 * kMoundRadiusMm / voxelSizeYMm);
    if (!std::isfinite(voxelSizeYMm)
        || !std::isfinite(requiredVoxelCountY)
        || requiredVoxelCountY > static_cast<double>(std::numeric_limits<int>::max() - 1)) {
        if (rank == 0) {
            std::cerr << "Error: Y voxel size produces an unsupported voxel count."
                      << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    int voxelCountY = std::max(2, static_cast<int>(requiredVoxelCountY));
    if (voxelCountY % 2 != 0) ++voxelCountY;
    const double moundBaseZMm = MoundBaseZMm(groundSlope);
    const double moundHeightMm = MoundHeightMm(groundSlope);
    const double voxelZMinMm =
        std::floor(moundBaseZMm / kVoxelSizeZMm) * kVoxelSizeZMm;
    const double voxelZMaxMm =
        std::ceil((moundBaseZMm + moundHeightMm) / kVoxelSizeZMm)
        * kVoxelSizeZMm;
    const int voxelCountZ =
        std::max(1, static_cast<int>(std::llround((voxelZMaxMm - voxelZMinMm)
                                                  / kVoxelSizeZMm)));

    if (rank == 0) {
        std::cout << "\n--- MPI Standalone EMML Tomography Reconstruction (No TV) ---\n" << std::endl;
        std::cout << "Reading from: " << inputFile << std::endl;
        if (skimFactor > 1) {
            std::cout << "SKIMMING ENABLED: Processing 1 out of every " << skimFactor << " events." << std::endl;
        }
        if (closureMode) {
            std::cout << "CLOSURE MODE: synthetic 10x central-cavity contrast, seed "
                      << closureSeed << std::endl;
        }
        if (diagnosticsOnly) {
            std::cout << "DIAGNOSTICS-ONLY MODE: no EMML update or output will be written."
                      << std::endl;
            std::cout << "Diagnostic reference lambda: " << diagnosticsLambda
                      << " rad^2/cm" << std::endl;
        }
        if (angleOnly) {
            std::cout << "ANGLE-ONLY MODE: excluding transverse-position and covariance terms."
                      << std::endl;
        }
        if (virtualFourPatch) {
            std::cout << "VIRTUAL FOUR-PATCH MODE: 4 m x 4 m X/Y detector pairs "
                      << "at 28.5 m radius; X pair z=4 m, Y pair z=9 m (-Y) to 3 m (+Y)."
                      << std::endl;
        }
        if (sparseFourPair) {
            std::cout << "SPARSE FOUR-PAIR MODE: four 6 m x 4 m opposing pairs "
                      << "at 0, 45, 90, and 135 degrees; global -Y-to-+Y "
                      << "hillside rise of 6 m; "
                      << "2 m voxels (27 x 27 x 7)." << std::endl;
        }
        if (sparseEightPair) {
            std::cout << "SPARSE EIGHT-PAIR MODE: eight 6 m x 4 m opposing pairs "
                      << "at 22.5-degree spacing; global -Y-to-+Y hillside "
                      << "rise of 6 m; "
                      << "2 m voxels (27 x 27 x 7); 192 detector units." << std::endl;
        }
        if (genericDetectorMode) {
            const double equivalentUnits = detectorPairs * detectorWidthM * detectorHeightM;
            std::cout << "VIRTUAL DETECTOR ARRAY MODE: " << detectorPairs
                      << " opposing pairs, " << detectorWidthM << " m x " << detectorHeightM
                      << " m panels at explicit angles from the tilted-frame +Y axis [";
            for (std::size_t i = 0; i < detectorPairAnglesDeg.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << detectorPairAnglesDeg[i];
            }
            std::cout << "] deg (positive toward +X), axial centre "
                      << detectorElevationM << " m above the lower cap, cylinder inclination of "
                      << detectorInclineDeg << " degrees; "
                      << "entry and exit may use any two distinct detector stacks; "
                      << "1 m X/Z voxels and " << voxelSizeYM
                      << " m Y voxels with ground-adaptive Z bounds; equivalent 2 m x 2 m units="
                      << equivalentUnits << "." << std::endl;
        }
        std::cout << "Maximum post-acceptance muons: "
                  << (maxAcceptedMuons > 0 ? std::to_string(maxAcceptedMuons) : "unlimited")
                  << std::endl;
        std::cout << "EMML initial lambda: " << initialLambda << " rad^2/cm" << std::endl;
        std::cout << "Ground inclination: " << groundInclineDeg
                  << " deg (+Y downhill, matching DetectorConstruction.cc)" << std::endl;
        std::cout << "Incoming horizontal direction cut: dy/sqrt(dx^2+dy^2) >= "
                  << minMuonHorizontalY
                  << (minMuonHorizontalY <= -1.0 ? " (disabled)" : "")
                  << std::endl;
        std::cout << "Voxel grid: 54 x " << voxelCountY << " x " << voxelCountZ
                  << ", voxel size X/Y/Z = 1 / " << voxelSizeYM
                  << " / 1 m, global Y=["
                  << -(voxelCountY * voxelSizeYM) / 2.0 << ", "
                  << (voxelCountY * voxelSizeYM) / 2.0 << "] m, global Z=["
                  << voxelZMinMm / 1000.0
                  << ", " << voxelZMaxMm / 1000.0 << "] m" << std::endl;
    }


    MoundTomographyEMML_MPI emml("EMML_TomographyResults");

    emml.SetGeometryParameters(54, voxelCountY, voxelCountZ,
                               kVoxelSizeXMm, voxelSizeYMm, kVoxelSizeZMm,
                               voxelZMinMm);
    emml.SetMoundParameters(kMoundRadiusMm, moundHeightMm, groundSlope,
                            moundBaseZMm);
    emml.SetAlgorithmParameters(20, initialLambda, 100.0);
    emml.SetAngleOnly(angleOnly);
    emml.SetMaxAcceptedMuons(maxAcceptedMuons);

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

    long long localXViewAccepted = 0;
    long long localYViewAccepted = 0;
    const int detectorPairCount = genericDetectorMode ? detectorPairs
                                : sparseEightPair ? 8 : (sparseFourPair ? 4 : 0);
    const double detectorWidthMm = genericDetectorMode ? detectorWidthM * 1000.0 : 6000.0;
    const double detectorHeightMm = genericDetectorMode ? detectorHeightM * 1000.0 : 4000.0;
    const int detectorPanelCount = 2 * detectorPairCount;
    std::vector<long long> localEntryPanelAccepted(detectorPanelCount, 0);
    std::vector<long long> localExitPanelAccepted(detectorPanelCount, 0);
    std::vector<char> entryPanelMatches(detectorPanelCount, 0);
    std::vector<char> exitPanelMatches(detectorPanelCount, 0);
    long long localCrossPanelAccepted = 0;
    long long localDirectionAccepted = 0;
    long long localDirectionRejected = 0;

    for (int i = startEvent; i < endEvent; ++i) {
        if (i % skimFactor != 0) continue;

        tree->GetEntry(i);

        TVector3 pIn1(in1_x, in1_y, in1_z);
        TVector3 pIn2(in2_x, in2_y, in2_z);
        TVector3 pOut1(out1_x, out1_y, out1_z);
        TVector3 pOut2(out2_x, out2_y, out2_z);

        // Select on incoming travel azimuth, independently of where the
        // uniformly distributed virtual detector panels are placed.
        const TVector3 incomingDelta = pIn2 - pIn1;
        const double incomingHorizontalMagnitude =
            std::hypot(incomingDelta.X(), incomingDelta.Y());
        if (incomingHorizontalMagnitude <= 1.0e-12
            || incomingDelta.Y() / incomingHorizontalMagnitude < minMuonHorizontalY) {
            ++localDirectionRejected;
            continue;
        }
        ++localDirectionAccepted;

        if (detectorPairCount > 0) {
            // Each opposing pair contributes two two-layer detector stacks.
            // A track may enter through any stack and leave through any other
            // stack: it is not restricted to the geometrically opposing mate.
            const double kHalfWidthMm = detectorWidthMm / 2.0;
            const double kHalfHeightMm = detectorHeightMm / 2.0;
            constexpr double kPi = 3.14159265358979323846;
            const double detectorAxialCentreMm = detectorElevationM * 1000.0;
            const double detectorIncline =
                detectorInclineDeg * kPi / 180.0;
            const double detectorCos = std::cos(detectorIncline);
            const double detectorSin = std::sin(detectorIncline);
            const auto localX = [](const TVector3& p) {
                return p.X();
            };
            const auto localY = [detectorCos, detectorSin](const TVector3& p) {
                return detectorCos*p.Y() - detectorSin*p.Z();
            };
            const auto localAxial =
                [detectorCos, detectorSin](const TVector3& p) {
                    return detectorSin*p.Y() + detectorCos*p.Z();
                };

            std::fill(entryPanelMatches.begin(), entryPanelMatches.end(), 0);
            std::fill(exitPanelMatches.begin(), exitPanelMatches.end(), 0);
            for (int pair = 0; pair < detectorPairCount; ++pair) {
                double c = 0.0;
                double s = 0.0;
                if (genericDetectorMode) {
                    // User azimuth is measured from tilted-frame +Y, with
                    // positive angles rotating toward +X. The radial unit
                    // vector is therefore (sin(angle), cos(angle)).
                    const double angle =
                        detectorPairAnglesDeg[pair] * kPi / 180.0;
                    c = std::sin(angle);
                    s = std::cos(angle);
                } else {
                    // Preserve the historical evenly spaced geometry for
                    // the named sparse-four/eight legacy presets.
                    const double phi = pair * kPi / detectorPairCount;
                    c = std::cos(phi);
                    s = std::sin(phi);
                }
                const auto u = [c, s, &localX, &localY](const TVector3& p) {
                    return localX(p) * c + localY(p) * s;
                };
                const auto v = [c, s, &localX, &localY](const TVector3& p) {
                    return -localX(p) * s + localY(p) * c;
                };
                for (int side = 0; side < 2; ++side) {
                    const double sign = side == 0 ? -1.0 : 1.0;
                    const int panel = 2 * pair + side;
                    const auto onPanel = [&](const TVector3& p) {
                        return sign * u(p) > 0.0 && std::abs(v(p)) <= kHalfWidthMm
                            && std::abs(localAxial(p) - detectorAxialCentreMm)
                               <= kHalfHeightMm;
                    };
                    entryPanelMatches[panel] = onPanel(pIn1) && onPanel(pIn2);
                    exitPanelMatches[panel] = onPanel(pOut1) && onPanel(pOut2);
                }
            }

            bool accepted = false;
            for (int entryPanel = 0; entryPanel < detectorPanelCount; ++entryPanel) {
                if (!entryPanelMatches[entryPanel]) continue;
                for (int exitPanel = 0; exitPanel < detectorPanelCount; ++exitPanel) {
                    if (entryPanel != exitPanel && exitPanelMatches[exitPanel]) {
                        accepted = true;
                        break;
                    }
                }
                if (accepted) break;
            }
            if (!accepted) continue;
            ++localCrossPanelAccepted;
            for (int panel = 0; panel < detectorPanelCount; ++panel) {
                if (entryPanelMatches[panel]) ++localEntryPanelAccepted[panel];
                if (exitPanelMatches[panel]) ++localExitPanelAccepted[panel];
            }
        } else if (virtualFourPatch) {
            // The ROOT hits come from full cylindrical layers at 28 and 29 m.
            // Emulate two opposing 4 m x 4 m detector patches in each view by
            // requiring both layers on the entry and exit sides to be crossed.
            constexpr double kHalfPatchMm = 2000.0;
            constexpr double kXViewZMm = 4000.0;
            constexpr double kYInZMm = 9000.0;   // -Y incoming detector (6 m above +Y)
            constexpr double kYOutZMm = 3000.0;  // +Y outgoing detector

            const auto inXPatch = [&](const TVector3& p) {
                return std::abs(p.Y()) <= kHalfPatchMm
                    && std::abs(p.Z() - kXViewZMm) <= kHalfPatchMm;
            };
            const auto inYPatch = [&](const TVector3& p, double zCenter) {
                return std::abs(p.X()) <= kHalfPatchMm
                    && std::abs(p.Z() - zCenter) <= kHalfPatchMm;
            };

            const bool xForward = pIn1.X() < 0.0 && pIn2.X() < 0.0
                               && pOut1.X() > 0.0 && pOut2.X() > 0.0
                               && inXPatch(pIn1) && inXPatch(pIn2)
                               && inXPatch(pOut1) && inXPatch(pOut2);
            const bool xReverse = pIn1.X() > 0.0 && pIn2.X() > 0.0
                               && pOut1.X() < 0.0 && pOut2.X() < 0.0
                               && inXPatch(pIn1) && inXPatch(pIn2)
                               && inXPatch(pOut1) && inXPatch(pOut2);
            const bool yForward = pIn1.Y() < 0.0 && pIn2.Y() < 0.0
                               && pOut1.Y() > 0.0 && pOut2.Y() > 0.0
                               && inYPatch(pIn1, kYInZMm) && inYPatch(pIn2, kYInZMm)
                               && inYPatch(pOut1, kYOutZMm) && inYPatch(pOut2, kYOutZMm);

            if (!xForward && !xReverse && !yForward) continue;
            if (xForward || xReverse) ++localXViewAccepted;
            if (yForward) ++localYViewAccepted;
        }

        MuonScatteringEvent ev;
        ev.EventId = i;
        
        ev.EntryPoint = pIn2; 
        ev.ExitPoint = pOut1; 
        
        ev.EntryDir = (pIn2 - pIn1).Unit(); 
        ev.ExitDir = (pOut2 - pOut1).Unit();  

	    
        // Project hits to the analytic ellipsoid clipped by the ground plane.
        if (!IntersectMound(pIn2, ev.EntryDir, groundSlope, true, ev.EntryPoint)
            || !IntersectMound(pOut1, ev.ExitDir, groundSlope, false, ev.ExitPoint)) {
            continue;
        }
 
        double p_magnitude = std::sqrt((in2_px * in2_px) + (in2_py * in2_py) + (in2_pz * in2_pz));
	ev.p = p_magnitude;
        
        ev.ScatteringAngle = ev.EntryDir.Angle(ev.ExitDir); 
        
        emml.AddEvent(ev);
    }

    long long globalDirectionAccepted = 0;
    long long globalDirectionRejected = 0;
    MPI_Allreduce(&localDirectionAccepted, &globalDirectionAccepted, 1,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localDirectionRejected, &globalDirectionRejected, 1,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Incoming direction selection: accepted=" << globalDirectionAccepted
                  << ", rejected=" << globalDirectionRejected << std::endl;
    }

    if (virtualFourPatch) {
        long long globalXViewAccepted = 0;
        long long globalYViewAccepted = 0;
        MPI_Allreduce(&localXViewAccepted, &globalXViewAccepted, 1, MPI_LONG_LONG,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localYViewAccepted, &globalYViewAccepted, 1, MPI_LONG_LONG,
                      MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Virtual detector acceptance: X view=" << globalXViewAccepted
                      << ", Y view=" << globalYViewAccepted
                      << ", total=" << globalXViewAccepted + globalYViewAccepted
                      << std::endl;
        }
    }
    if (detectorPairCount > 0) {
        long long globalCrossPanelAccepted = 0;
        std::vector<long long> globalEntryPanelAccepted(detectorPanelCount, 0);
        std::vector<long long> globalExitPanelAccepted(detectorPanelCount, 0);
        MPI_Allreduce(&localCrossPanelAccepted, &globalCrossPanelAccepted, 1,
                      MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(localEntryPanelAccepted.data(), globalEntryPanelAccepted.data(),
                      detectorPanelCount, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(localExitPanelAccepted.data(), globalExitPanelAccepted.data(),
                      detectorPanelCount, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Virtual detector array acceptance: unique cross-stack events="
                      << globalCrossPanelAccepted << ", entry stacks=";
            for (int panel = 0; panel < detectorPanelCount; ++panel) {
                if (panel > 0) std::cout << ", ";
                std::cout << globalEntryPanelAccepted[panel];
            }
            std::cout << ", exit stacks=";
            for (int panel = 0; panel < detectorPanelCount; ++panel) {
                if (panel > 0) std::cout << ", ";
                std::cout << globalExitPanelAccepted[panel];
            }
            std::cout << std::endl;
        }
    }
    
    fin->Close();

    // Run processing
    emml.TracePaths(); 
    if (closureMode) {
        emml.GenerateCentralCavityClosure(1.0e-7, 0.1, closureSeed);
    }
    if (diagnosticsOnly) {
        emml.PrintResidualDiagnostics(diagnosticsLambda);
        MPI_Finalize();
        return 0;
    }
    emml.RunEMML();
    
    if (rank == 0) {
        // 1. Save EVERYTHING to the ROOT file first
        emml.SaveResults(-20.0, 0.0, 1); 
    }
    MPI_Finalize();
    return 0;
}
