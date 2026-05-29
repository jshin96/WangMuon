# WangMuon: Cosmic Ray Geant4 Simulation

## Overview
WangMuon is a Geant4-based Monte Carlo simulation designed to model high-energy atmospheric cosmic muons interacting with geological structures. 

Specifically, it simulates 10–100 GeV muons raining down on a 10-meter high, 20-meter radius dirt mound ($\text{SiO}_2$) that contains a hidden 3x3x3 meter empty room at its center. The simulation tracks energy deposition, boundary crossings, and generates both an event-by-event ROOT ntuple for analysis and automatic 3D track visualizations.

## Project Structure & Source Code
Here is a brief explanation of what each core class does in the simulation:

* **`main.cc`**: The entry point of the application. It establishes the `G4RunManager` (forced to sequential mode to allow GUI screenshot commands) and initializes the Qt visualizer.
* **`DetectorConstruction.cc`**: Builds the physical geometry. It constructs the invisible World box, the semi-transparent Dirt Mound (`G4Ellipsoid`), the solid Ground plane, and the blue Empty Room (`G4Box`). Applies custom `G4VisAttributes` to each.
* **`PrimaryGeneratorAction.cc`**: The cosmic ray gun. Generates primary muons with random energies (10-100 GeV). Positions them uniformly on a 20m disk above the mound and fires them downward using a realistic $\cos(\theta)$ zenith angular distribution and a uniform azimuthal distribution.
* **`RunAction.cc`**: Sets up the ROOT analysis manager. Defines the structure of `wang_muon_data.root`, creating scalar columns for entry/exit coordinates and momenta, and vector columns for step-by-step tracking. 
* **`EventAction.cc`**: Resets the data variables at the beginning of each event. At the end of the event, it commits the scalar data to the ROOT file and uses the `G4UImanager` to automatically capture and rename a `.jpg` screenshot of the 3D event track.
* **`SteppingAction.cc`**: The active tracker. It acts as a boundary detective, recording the exact 4-momentum and position when the primary muon enters and exits the mound/ground. It also flags if the muon passed through the room and records its position/momentum anytime its kinetic energy changes.
* **`TrackingAction.cc`**: Currently bypassed (logic handled in SteppingAction).
* **`ActionInitialization.cc`**: The registry that binds all the user actions (Generator, Run, Event, Stepping) to the Geant4 Run Manager.

## Prerequisites
To compile and run this project, you need:
* **Geant4** (Compiled with Qt5/Qt6 and OpenGL support)
* **ROOT** (For `wang_muon_data.root` ntuple analysis)
* **CMake**
* **C++ Compiler** (Clang/GCC)
```bash
# 1. Create a new environment named WangMuon and install the required packages
conda create -n WangMuon -c conda-forge geant4 root cmake compilers

# 2. Activate the environment
conda activate WangMuon


## How to Compile
This project uses standard CMake out-of-source building. Open your terminal in the main project directory and run:

```bash
# 1. Create a build directory
mkdir build
cd build

# 2. Generate the Makefiles
cmake ..

# 3. Compile the code
make

#Once the GUI opens, type 

/run/beamOn 10

#to generate 10 events.
