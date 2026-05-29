# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(WangMuon)

# Find the Geant4 package
find_package(Geant4 REQUIRED ui_all vis_all)

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

# Tell CMake to compile our main.cc file into an executable named "wang_muon"
add_executable(wang_muon 
    main.cc 
    DetectorConstruction.cc 
    PrimaryGeneratorAction.cc 
    ActionInitialization.cc
    EventAction.cc
    SteppingAction.cc
    RunAction.cc
    TrackingAction.cc
)

# Link the Geant4 libraries to our executable
target_link_libraries(wang_muon ${Geant4_LIBRARIES})
