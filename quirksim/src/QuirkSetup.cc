#include "QuirkSetup.hh"
#include "QSimConfig.hh"

#include "G4Exception.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#ifdef QSIM_WITH_PYTHIA8
#include "Pythia8/Pythia.h"
#endif

QuirkConfig& quirkConfig() {
  static QuirkConfig cfg;
  return cfg;
}

std::vector<std::string> quirkPythiaMassSettings(double m) {
  // the pair needs at least 2 m; Pythia DY of the lepton through the Z' machinery
  return {"17:m0 = " + std::to_string(m), "17:mWidth = 0.", "17:mayDecay = off", "17:isResonance = false",
          "32:mMin = " + std::to_string(2. * m + 1.)};
}

#ifdef QSIM_WITH_PYTHIA8
std::unique_ptr<Pythia8::Pythia> makeQuirkPythia(int seedOffset) {
  const QSimConfig& c = qsimConfig();
  auto p = std::make_unique<Pythia8::Pythia>("../share/Pythia8/xmldoc", false);  // PYTHIA8DATA wins
  if (!c.pythiaCmnd.empty() && !p->readFile(c.pythiaCmnd))
    G4Exception("makeQuirkPythia", "QSimPythia", FatalException, ("cannot read " + c.pythiaCmnd).c_str());
  for (const auto& s : c.pythiaStrings)
    p->readString(s);
  for (const auto& s : quirkPythiaMassSettings(quirkConfig().mass / CLHEP::GeV))
    p->readString(s);
  p->readString("Random:setSeed = on");
  p->readString("Random:seed = " + std::to_string(c.pythiaSeed + seedOffset));
  p->readString("Next:numberCount = 0");
  if (!p->init())
    G4Exception("makeQuirkPythia", "QSimPythia", FatalException, "Pythia8 init failed");
  return p;
}
#endif
