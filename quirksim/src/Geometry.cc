#include "Geometry.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4AutoDelete.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <map>
#include <sstream>

GeomInfo& geomInfo() {
  static GeomInfo g;
  return g;
}

G4int GeomInfo::detIndex(const std::string& name) {
  for (size_t i = 0; i < dets.size(); ++i)
    if (dets[i] == name)
      return i;
  dets.push_back(name);
  return dets.size() - 1;
}

void QSimField::GetFieldValue(const G4double point[4], G4double* b) const {
  b[0] = b[1] = b[2] = 0.;
  const G4double r = std::hypot(point[0], point[1]);
  const G4double az = std::abs(point[2]);
  for (const auto& f : m_regions) {
    if (r >= f.rmin && r < f.rmax && az >= f.zmin && az < f.zmax) {
      b[2] = f.bz;
      return;
    }
  }
}

namespace {
  [[noreturn]] void cardError(const std::string& card, int line, const std::string& msg) {
    G4ExceptionDescription ed;
    ed << card << ":" << line << ": " << msg;
    G4Exception("DetectorConstruction", "QSimCard", FatalException, ed);
    throw;  // not reached
  }
}  // namespace

G4Material* DetectorConstruction::material(const std::string& name) {
  if (auto m = G4Material::GetMaterial(name, false))
    return m;
  return G4NistManager::Instance()->FindOrBuildMaterial(name);
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  GeomInfo& g = geomInfo();
  g = GeomInfo();
  std::ifstream in(m_card);
  if (!in)
    cardError(m_card, 0, "cannot open card");

  G4LogicalVolume* worldLV = nullptr;
  G4VPhysicalVolume* worldPV = nullptr;
  G4bool checkOverlaps = true;
  std::string line;
  int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    if (auto h = line.find('#'); h != std::string::npos)
      line.erase(h);
    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key))
      continue;
    if (key == "world") {
      std::string mat;
      ss >> g.worldR >> g.worldZ >> mat;
      g.worldR *= cm;
      g.worldZ *= cm;
      auto m = material(mat);
      if (!m)
        cardError(m_card, ln, "unknown material " + mat);
      worldLV = new G4LogicalVolume(new G4Tubs("world", 0, g.worldR, g.worldZ, 0, CLHEP::twopi), m, "world");
      worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0);
      VolInfo w;
      w.name = "world";
      w.det = g.detIndex("WORLD");
      w.rmax = g.worldR;
      w.zmin = -g.worldZ;
      w.zmax = g.worldZ;
      g.vols.push_back(w);
    } else if (key == "checkOverlaps") {
      ss >> checkOverlaps;
    } else if (key == "mix") {
      // mix <name> <density g/cm3> <nist> <mass fraction> ...
      std::string name, comp;
      double rho, frac;
      ss >> name >> rho;
      std::vector<std::pair<G4Material*, double>> parts;
      double sum = 0;
      while (ss >> comp >> frac) {
        auto m = material(comp);
        if (!m)
          cardError(m_card, ln, "unknown material " + comp);
        parts.emplace_back(m, frac);
        sum += frac;
      }
      auto m = new G4Material(name, rho * CLHEP::g / cm3, parts.size());
      for (auto& [pm, f] : parts)
        m->AddMaterial(pm, f / sum);
    } else if (key == "scaled") {
      // scaled <name> <nist> <density factor>
      std::string name, base;
      double f;
      ss >> name >> base >> f;
      auto bm = material(base);
      if (!bm)
        cardError(m_card, ln, "unknown material " + base);
      G4NistManager::Instance()->BuildMaterialWithNewDensity(name, base, bm->GetDensity() * f);
    } else if (key == "tube") {
      // tube <name> <mat> <rmin> <rmax|+x0:v|+lam:v|+cm:v> <zmin> <zmax|+...> [det=] [layer=] [sens] [calo] [mirror]
      std::string name, mat, srmax, szmax, opt;
      double rmin, zmin;
      ss >> name >> mat >> rmin >> srmax >> zmin >> szmax;
      if (!worldLV)
        cardError(m_card, ln, "world must come first");
      auto m = material(mat);
      if (!m)
        cardError(m_card, ln, "unknown material " + mat);
      auto extent = [&](double lo, const std::string& s) {
        if (s.empty())
          cardError(m_card, ln, "missing extent");
        if (s[0] != '+')
          return std::stod(s) * cm;
        auto c = s.find(':');
        const std::string unit = s.substr(1, c - 1);
        const double v = std::stod(s.substr(c + 1));
        if (unit == "x0")
          return lo * cm + v * m->GetRadlen();
        if (unit == "lam")
          return lo * cm + v * m->GetNuclearInterLength();
        if (unit == "cm")
          return lo * cm + v * cm;
        if (unit == "um")
          return lo * cm + v * um;
        cardError(m_card, ln, "unknown thickness unit " + unit);
      };
      VolInfo v;
      v.name = name;
      v.rmin = rmin * cm;
      v.rmax = extent(rmin, srmax);
      v.zmin = zmin * cm;
      v.zmax = extent(zmin, szmax);
      bool mirror = false;
      std::string det = name;
      while (ss >> opt) {
        if (opt.rfind("det=", 0) == 0)
          det = opt.substr(4);
        else if (opt.rfind("layer=", 0) == 0)
          v.layer = std::stoi(opt.substr(6));
        else if (opt == "sens")
          v.sens = true;
        else if (opt == "calo")
          v.calo = true;
        else if (opt == "mirror")
          mirror = true;
        else
          cardError(m_card, ln, "unknown option " + opt);
      }
      if (v.rmax <= v.rmin || v.zmax <= v.zmin)
        cardError(m_card, ln, "empty volume " + name);
      v.det = g.detIndex(det);
      auto lv = new G4LogicalVolume(
          new G4Tubs(name, v.rmin, v.rmax, 0.5 * (v.zmax - v.zmin), 0, CLHEP::twopi), m, name);
      for (int side = 0; side < (mirror ? 2 : 1); ++side) {
        VolInfo vi = v;
        if (side == 1) {
          vi.zmin = -v.zmax;
          vi.zmax = -v.zmin;
          vi.name += "_m";
        }
        const G4int copy = g.vols.size();
        g.vols.push_back(vi);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0.5 * (vi.zmin + vi.zmax)), lv, vi.name, worldLV, false, copy,
                          checkOverlaps);
      }
    } else if (key == "field") {
      FieldRegion f;
      ss >> f.bz >> f.rmin >> f.rmax >> f.zmin >> f.zmax;
      f.bz *= tesla;
      f.rmin *= cm;
      f.rmax *= cm;
      f.zmin *= cm;
      f.zmax *= cm;
      g.field.push_back(f);
    } else if (key == "boundary") {
      Boundary b;
      ss >> b.name >> b.r >> b.halfz;
      b.r *= cm;
      b.halfz *= cm;
      g.boundaries.push_back(b);
    } else {
      cardError(m_card, ln, "unknown keyword " + key);
    }
  }
  if (!worldPV)
    cardError(m_card, ln, "no world");
  G4cout << "DetectorConstruction: " << m_card << ": " << g.vols.size() - 1 << " volumes, " << g.dets.size()
         << " detector groups, " << g.field.size() << " field regions" << G4endl;
  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  // one field and field manager per thread
  const GeomInfo& g = geomInfo();
  if (g.field.empty())
    return;
  auto field = new QSimField(g.field);
  G4AutoDelete::Register(field);
  auto fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(field);
  fieldMgr->CreateChordFinder(field);
}
