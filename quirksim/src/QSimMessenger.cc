#include "QSimMessenger.hh"
#include "QSimConfig.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

QSimConfig& qsimConfig() {
  static QSimConfig cfg;
  return cfg;
}

namespace {
  template <class T>
  T* setup(T* c, const char* guide, bool preInitOnly = false) {
    c->SetGuidance(guide);
    if (preInitOnly)
      c->AvailableForStates(G4State_PreInit);
    else
      c->AvailableForStates(G4State_PreInit, G4State_Idle);
    c->SetToBeBroadcasted(false);
    return c;
  }
}  // namespace

QSimMessenger::QSimMessenger() {
  for (const char* d : {"/qsim/", "/qsim/gen/", "/qsim/gun/", "/qsim/pythia/", "/qsim/out/"}) {
    m_dirs.emplace_back(new G4UIdirectory(d));
  }
  auto str = [this](const char* p, const char* g, bool pre = false) {
    m_cmds.emplace_back(setup(new G4UIcmdWithAString(p, this), g, pre));
  };
  auto dim = [this](const char* p, const char* g, const char* unit, bool pre = false) {
    auto c = setup(new G4UIcmdWithADoubleAndUnit(p, this), g, pre);
    c->SetDefaultUnit(unit);
    m_cmds.emplace_back(c);
  };
  auto num = [this](const char* p, const char* g) { m_cmds.emplace_back(setup(new G4UIcmdWithADouble(p, this), g)); };
  auto integer = [this](const char* p, const char* g) {
    m_cmds.emplace_back(setup(new G4UIcmdWithAnInteger(p, this), g));
  };
  auto flag = [this](const char* p, const char* g) { m_cmds.emplace_back(setup(new G4UIcmdWithABool(p, this), g)); };
  auto vec = [this](const char* p, const char* g, const char* unit) {
    auto c = setup(new G4UIcmdWith3VectorAndUnit(p, this), g);
    c->SetDefaultUnit(unit);
    m_cmds.emplace_back(c);
  };

  str("/qsim/card", "Detector card file.", true);
  str("/qsim/physList", "Reference physics list name (G4PhysListFactory).", true);
  dim("/qsim/prodCut", "Default production cut.", "mm", true);
  dim("/qsim/maxTime", "Kill all tracks after this global time (<= 0: never).", "ns");
  flag("/qsim/killSecondaries", "Kill all secondaries (quirks only).");

  str("/qsim/gen/mode", "gun | hepmc | pythia");
  str("/qsim/gen/hepmcFile", "HepMC3 ascii input file.");
  flag("/qsim/gen/quirksOnly", "Pass only the quirks to Geant4.");
  num("/qsim/gen/etaMax", "Skip generator particles beyond this |eta|.");
  integer("/qsim/gen/skipEvents", "Skip events at the start of the HepMC3 file.");
  dim("/qsim/gen/vtxSigmaXY", "Gaussian vertex smearing in x and y.", "mm");
  dim("/qsim/gen/vtxSigmaZ", "Gaussian vertex smearing in z.", "mm");

  vec("/qsim/gun/p1", "Momentum of the quirk (positive PDG code).", "GeV");
  vec("/qsim/gun/p2", "Momentum of the antiquirk.", "GeV");
  vec("/qsim/gun/vertex", "Production point.", "mm");

  str("/qsim/pythia/cmnd", "Pythia8 settings file.");
  str("/qsim/pythia/readString", "One Pythia8 setting.");
  integer("/qsim/pythia/seed", "Pythia8 seed (plus the thread id).");

  str("/qsim/out/file", "Output ROOT file.");
  integer("/qsim/out/stepPrescale", "Write every N-th quirk step (0: none).");
  flag("/qsim/out/towers", "Write calorimeter towers.");
  integer("/qsim/out/printEvery", "Print a line every N events.");
}

QSimMessenger::~QSimMessenger() = default;

void QSimMessenger::SetNewValue(G4UIcommand* cmd, G4String v) {
  QSimConfig& c = qsimConfig();
  const G4String& p = cmd->GetCommandPath();
  auto d = [&]() { return G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(v); };
  auto b = [&]() { return G4UIcommand::ConvertToBool(v); };
  auto i = [&]() { return G4UIcommand::ConvertToInt(v); };
  auto v3 = [&]() { return G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(v); };
  if (p == "/qsim/card")
    c.card = v;
  else if (p == "/qsim/physList")
    c.physList = v;
  else if (p == "/qsim/prodCut")
    c.prodCut = d();
  else if (p == "/qsim/maxTime")
    c.maxTime = d();
  else if (p == "/qsim/killSecondaries")
    c.killSecondaries = b();
  else if (p == "/qsim/gen/mode")
    c.genMode = v;
  else if (p == "/qsim/gen/hepmcFile")
    c.hepmcFile = v;
  else if (p == "/qsim/gen/quirksOnly")
    c.quirksOnly = b();
  else if (p == "/qsim/gen/etaMax")
    c.etaMax = G4UIcommand::ConvertToDouble(v);
  else if (p == "/qsim/gen/skipEvents")
    c.skipEvents = i();
  else if (p == "/qsim/gen/vtxSigmaXY")
    c.vtxSigmaXY = d();
  else if (p == "/qsim/gen/vtxSigmaZ")
    c.vtxSigmaZ = d();
  else if (p == "/qsim/gun/p1")
    c.gunP1 = v3();
  else if (p == "/qsim/gun/p2")
    c.gunP2 = v3();
  else if (p == "/qsim/gun/vertex")
    c.gunVertex = v3();
  else if (p == "/qsim/pythia/cmnd")
    c.pythiaCmnd = v;
  else if (p == "/qsim/pythia/readString")
    c.pythiaStrings.push_back(v);
  else if (p == "/qsim/pythia/seed")
    c.pythiaSeed = i();
  else if (p == "/qsim/out/file")
    c.output = v;
  else if (p == "/qsim/out/stepPrescale")
    c.stepPrescale = i();
  else if (p == "/qsim/out/towers")
    c.towers = b();
  else if (p == "/qsim/out/printEvery")
    c.printEvery = i();
}
