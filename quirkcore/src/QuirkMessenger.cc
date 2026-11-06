#include "quirkcore/QuirkMessenger.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

namespace {
  template <class T>
  T* withStates(T* c) {
    c->AvailableForStates(G4State_PreInit);
    c->SetToBeBroadcasted(false);
    return c;
  }
  G4UIcmdWithADoubleAndUnit* dimCmd(
      G4UImessenger* m, const char* path, const char* guide, const char* unit, const char* cat) {
    auto c = new G4UIcmdWithADoubleAndUnit(path, m);
    c->SetGuidance(guide);
    c->SetUnitCategory(cat);
    c->SetDefaultUnit(unit);
    return withStates(c);
  }
}  // namespace

QuirkMessenger::QuirkMessenger(QuirkConfig& cfg) : m_cfg(cfg) {
  m_dir = std::make_unique<G4UIdirectory>("/quirk/");
  m_dir->SetGuidance("Quirk parameters (defaults of g4SimHits_cfi.py); PreInit only.");

  auto add = [this](G4UIcommand* c) { m_cmds.emplace_back(c); };
  add(dimCmd(this, "/quirk/mass", "Quirk mass (> 0 enables quirks).", "GeV", "Energy"));
  add(dimCmd(this, "/quirk/lambda", "Infracolor scale Lambda, string tension Lambda^2/(hbar c).", "eV", "Energy"));
  add(dimCmd(this, "/quirk/firstStringLength", "Length of the first string vector.", "mm", "Length"));
  add(dimCmd(this, "/quirk/maxMerge", "Max merge length of string vectors.", "mm", "Length"));
  add(dimCmd(this, "/quirk/maxRadius", "Envelope radius, quirks killed outside (<= 0: never).", "cm", "Length"));
  add(dimCmd(this, "/quirk/maxZ", "Envelope half length, quirks killed outside (<= 0: never).", "cm", "Length"));
  add(dimCmd(this, "/quirk/stopThreshold", "Partner kinetic energy below which the pair counts as stopped.", "MeV", "Energy"));
  add(dimCmd(this, "/quirk/crossingLength", "Restart the string below this length (0: never, experimental).", "mm", "Length"));
  {
    auto c = withStates(new G4UIcmdWithADouble("/quirk/charge", this));
    c->SetGuidance("Electric charge of the positive PDG code, in units of e.");
    add(c);
  }
  {
    auto c = withStates(new G4UIcmdWithADouble("/quirk/maxBoost", this));
    c->SetGuidance("Max rapidity change per step (CMSSW 0.1, Athena 0.01).");
    c->SetParameterName("maxBoost", false);
    c->SetRange("maxBoost>0 && maxBoost<1");
    add(c);
  }
  {
    auto c = withStates(new G4UIcmdWithAnInteger("/quirk/pdg", this));
    c->SetGuidance("PDG code of the quirk.");
    add(c);
  }
  {
    auto c = withStates(new G4UIcmdWithAnInteger("/quirk/looperTrials", this));
    c->SetGuidance("Looper trials before killing.");
    add(c);
  }
  {
    auto c = withStates(new G4UIcmdWithAnInteger("/quirk/verbose", this));
    c->SetGuidance("Verbosity.");
    add(c);
  }
  for (const char* name : {"/quirk/keepStopped", "/quirk/moveTinySteps", "/quirk/pairProduction"}) {
    auto c = withStates(new G4UIcmdWithABool(name, this));
    c->SetGuidance("Boolean switch, see QuirkConfig.h.");
    add(c);
  }
}

QuirkMessenger::~QuirkMessenger() = default;

void QuirkMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (m_locked) {
    G4Exception("QuirkMessenger::SetNewValue", "QuirkLocked", JustWarning,
                ("quirk particles already built, ignoring " + cmd->GetCommandPath()).c_str());
    return;
  }
  const G4String& n = cmd->GetCommandName();
  auto d = [&]() { return G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(value); };
  if (n == "mass")
    m_cfg.mass = d();
  else if (n == "lambda")
    m_cfg.lambda = d();
  else if (n == "firstStringLength")
    m_cfg.firstStringLength = d();
  else if (n == "maxMerge")
    m_cfg.maxMerge = d();
  else if (n == "maxRadius")
    m_cfg.maxRadius = d();
  else if (n == "maxZ")
    m_cfg.maxZ = d();
  else if (n == "stopThreshold")
    m_cfg.stopThreshold = d();
  else if (n == "crossingLength")
    m_cfg.crossingLength = d();
  else if (n == "charge")
    m_cfg.charge = G4UIcommand::ConvertToDouble(value) * CLHEP::eplus;
  else if (n == "maxBoost")
    m_cfg.maxBoost = G4UIcommand::ConvertToDouble(value);
  else if (n == "pdg")
    m_cfg.pdg = G4UIcommand::ConvertToInt(value);
  else if (n == "looperTrials")
    m_cfg.looperTrials = G4UIcommand::ConvertToInt(value);
  else if (n == "verbose")
    m_cfg.verbose = G4UIcommand::ConvertToInt(value);
  else if (n == "keepStopped")
    m_cfg.keepStopped = G4UIcommand::ConvertToBool(value);
  else if (n == "moveTinySteps")
    m_cfg.moveTinySteps = G4UIcommand::ConvertToBool(value);
  else if (n == "pairProduction")
    m_cfg.pairProduction = G4UIcommand::ConvertToBool(value);
}

G4String QuirkMessenger::GetCurrentValue(G4UIcommand* cmd) {
  const G4String& n = cmd->GetCommandName();
  auto u = [&](G4double v, const char* unit) { return G4UIcommand::ConvertToString(v / G4UIcommand::ValueOf(unit)) + " " + unit; };
  if (n == "mass")
    return u(m_cfg.mass, "GeV");
  if (n == "lambda")
    return u(m_cfg.lambda, "eV");
  if (n == "firstStringLength")
    return u(m_cfg.firstStringLength, "mm");
  if (n == "maxMerge")
    return u(m_cfg.maxMerge, "mm");
  if (n == "maxRadius")
    return u(m_cfg.maxRadius, "cm");
  if (n == "maxZ")
    return u(m_cfg.maxZ, "cm");
  if (n == "stopThreshold")
    return u(m_cfg.stopThreshold, "MeV");
  if (n == "crossingLength")
    return u(m_cfg.crossingLength, "mm");
  if (n == "charge")
    return G4UIcommand::ConvertToString(m_cfg.charge / CLHEP::eplus);
  if (n == "maxBoost")
    return G4UIcommand::ConvertToString(m_cfg.maxBoost);
  if (n == "pdg")
    return G4UIcommand::ConvertToString(m_cfg.pdg);
  if (n == "looperTrials")
    return G4UIcommand::ConvertToString(m_cfg.looperTrials);
  if (n == "verbose")
    return G4UIcommand::ConvertToString(m_cfg.verbose);
  if (n == "keepStopped")
    return G4UIcommand::ConvertToString(m_cfg.keepStopped);
  if (n == "moveTinySteps")
    return G4UIcommand::ConvertToString(m_cfg.moveTinySteps);
  if (n == "pairProduction")
    return G4UIcommand::ConvertToString(m_cfg.pairProduction);
  return "";
}
