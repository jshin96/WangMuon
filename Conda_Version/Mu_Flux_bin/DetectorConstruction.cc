#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Exception.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4VisAttributes.hh"
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
namespace {
// Read a number from the job settings.  If there is no setting, use the
// fallback number written in the code instead.
G4double Env(const char* name, G4double fallback) {
  const char* text = std::getenv(name);
  if (!text) return fallback;

  try {
    std::size_t used = 0;
    const auto value = std::stod(text, &used);
    if (used != std::string(text).size() || !std::isfinite(value)) {
      throw std::invalid_argument("not a finite number");
    }
    return value;
  } catch (...) {
    G4Exception("DetectorConstruction", "InvalidConfiguration", FatalException,
                (std::string(name) + " must be finite.").c_str());
    return fallback;
  }
}
}  // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // Air is the empty space.  The GEM active gas is 70% argon and 30% CO2 by
  // volume.  Geant4 needs mass fractions, so 70/30 by volume becomes about
  // 68/32 by mass after accounting for the different molecule masses.
  auto* nist = G4NistManager::Instance();
  auto* air = nist->FindOrBuildMaterial("G4_AIR");
  auto* argon = nist->FindOrBuildMaterial("G4_Ar");
  auto* carbonDioxide = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  auto* gas = G4Material::GetMaterial("GEMGasArCO2", false);
  if (!gas) {
    gas = new G4Material("GEMGasArCO2", 1.84 * mg/cm3, 2,
                         kStateGas, 293.15 * kelvin, atmosphere);
    gas->AddMaterial(argon, 0.680);
    gas->AddMaterial(carbonDioxide, 0.320);
  }
  auto* world = new G4LogicalVolume(
      new G4Box("World", 10*m, 10*m, 10*m), air, "World");
  auto* pv = new G4PVPlacement(nullptr, {}, world, "World", nullptr,
                                false, 0, true);
  // IN2 is nearest the magnet; defaults are -50 and -70 cm (20 cm spacing).
  const auto in2 = Env("MUON_GEM_IN2_X_M", -0.50) * m;
  const auto in1 = Env("MUON_GEM_IN1_X_M", -0.70) * m;
  const auto out = Env("MUON_GEM_OUT_X_M", 1.00) * m;
  const auto size = Env("MUON_GEM_SIZE_CM", 10.0) * cm;
  const auto thick = Env("MUON_GEM_THICKNESS_MM", 3.0) * mm;
  if (!(in1 < in2 && in2 < 0.0 && out > 0.0 && size > 0.0 && thick > 0.0)) {
    G4Exception("DetectorConstruction", "InvalidGEMLayout", FatalException,
                "Require IN1_X < IN2_X < 0, OUT_X > 0, and positive GEM dimensions.");
  }
  // Make a thin, square detector board.  Its skinny side is X, so it faces -X.
  const auto place = [&](const char* name, G4double x) {
    auto* logical = new G4LogicalVolume(
        new G4Box(G4String(name) + "Solid", thick/2, size/2, size/2),
        gas, name);
    new G4PVPlacement(nullptr, {x, 0.0, 0.0}, logical, name, world,
                      false, 0, true);
    auto* attributes = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.6));
    attributes->SetForceSolid(true);
    logical->SetVisAttributes(attributes);
  };
  place("GEMIn1", in1);
  place("GEMIn2", in2);
  place("GEMOut", out);
  // Helmholtz-pair central region: uniform 1.5 T field along global +Y.
  const auto fieldT = Env("MUON_HELMHOLTZ_FIELD_T", 1.5);
  const auto length = Env("MUON_HELMHOLTZ_LENGTH_M", 0.70) * m;
  const auto radius = Env("MUON_HELMHOLTZ_APERTURE_CM", 30.0) * cm;
  if (fieldT < 0.0 || length <= 0.0 || radius <= 0.0) {
    G4Exception("DetectorConstruction", "InvalidField", FatalException,
                "Field magnitude must be nonnegative; length and radius positive.");
  }
  // This box is the useful middle of a Helmholtz pair.  Inside it the field
  // is uniform; outside it the field is zero in this simplified model.
  auto* fl = new G4LogicalVolume(
      new G4Tubs("HelmholtzFieldSolid", 0.0, radius, length/2,0, 2*CLHEP::pi),
//      new G4Box("HelmholtzFieldSolid", length/2, radius, radius),
      air, "HelmholtzField");
  new G4PVPlacement(nullptr, {}, fl, "HelmholtzField", world, false, 0, true);
  // (0, +B, 0) means the magnetic arrow points along the +Y direction.
  auto* field = new G4UniformMagField(G4ThreeVector{0.0, fieldT*tesla, 0.0});
  auto* fieldManager = new G4FieldManager(field);
  fieldManager->CreateChordFinder(field);
  fl->SetFieldManager(fieldManager, true);
  auto* fieldAttributes = new G4VisAttributes(G4Colour(1.0, 0.6, 0.0, 0.15));
  fieldAttributes->SetForceWireframe(true);
  fl->SetVisAttributes(fieldAttributes);
  world->SetVisAttributes(G4VisAttributes::GetInvisible());
  return pv;
}
