#ifndef QSIM_Geometry_hh
#define QSIM_Geometry_hh

// Parametric cylindrical detector from a text card (see cards/README.md),
// and the volume table shared read-only by all threads.

#include "G4VUserDetectorConstruction.hh"
#include "G4MagneticField.hh"
#include "globals.hh"

#include <string>
#include <vector>

class G4Material;

struct VolInfo {
  std::string name;
  G4int det{-1};  // index into GeomInfo::dets
  G4int layer{0};
  G4bool sens{false};  // quirk hits are recorded
  G4bool calo{false};  // energy goes to towers
  G4double rmin{0}, rmax{0}, zmin{0}, zmax{0};
};

struct FieldRegion {
  G4double bz, rmin, rmax, zmin, zmax;  // |z| range
};

struct Boundary {
  std::string name;
  G4double r, halfz;
};

struct GeomInfo {
  std::vector<VolInfo> vols;  // index = copy number of the placement, 0 = world
  std::vector<std::string> dets;
  std::vector<FieldRegion> field;
  std::vector<Boundary> boundaries;
  G4double worldR{0}, worldZ{0};
  G4int detIndex(const std::string& name);
};

GeomInfo& geomInfo();

class QSimField : public G4MagneticField {
public:
  explicit QSimField(const std::vector<FieldRegion>& regions) : m_regions(regions) {}
  void GetFieldValue(const G4double point[4], G4double* bfield) const override;

private:
  std::vector<FieldRegion> m_regions;
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction(const std::string& card) : m_card(card) {}
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

private:
  G4Material* material(const std::string& name);
  std::string m_card;
};

#endif
