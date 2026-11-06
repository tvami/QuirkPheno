// quirksim: quirk pairs through a parametric cylindrical detector (Geant4, no CMSSW).
//
// usage: quirksim -m run.mac [-t threads] [-n events] [-s seed] [-o out.root] [--mt]
//
// The macro lines before /run/initialize configure the job (/quirk/, /qsim/,
// /run/numberOfThreads ...); the particle definitions need the quirk mass, so
// the physics list is built only after them.

#include "Actions.hh"
#include "Geometry.hh"
#include "QSimConfig.hh"
#include "QSimMessenger.hh"
#include "QuirkSetup.hh"

#include "quirkcore/QuirkMessenger.h"
#include "quirkcore/QuirkPhysics.h"

#include "G4MTRunManager.hh"
#include "G4TaskRunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4UImanager.hh"
#include "G4Version.hh"
#include "Randomize.hh"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {
  void usage(const char* a) {
    G4cout << "usage: " << a << " -m macro [-t threads] [-n events] [-s seed] [-o out.root] [--mt]" << G4endl;
  }
}  // namespace

int main(int argc, char** argv) {
  std::string macro, output;
  int nthreads = 1, nevents = -1;
  long seed = 0;
  bool mt = false;
  for (int i = 1; i < argc; ++i) {
    auto next = [&]() -> const char* {
      if (i + 1 >= argc) {
        usage(argv[0]);
        std::exit(1);
      }
      return argv[++i];
    };
    if (!std::strcmp(argv[i], "-m"))
      macro = next();
    else if (!std::strcmp(argv[i], "-t"))
      nthreads = std::atoi(next());
    else if (!std::strcmp(argv[i], "-n"))
      nevents = std::atoi(next());
    else if (!std::strcmp(argv[i], "-s"))
      seed = std::atol(next());
    else if (!std::strcmp(argv[i], "-o"))
      output = next();
    else if (!std::strcmp(argv[i], "--mt"))
      mt = true;
    else {
      usage(argv[0]);
      return 1;
    }
  }
  if (macro.empty()) {
    usage(argv[0]);
    return 1;
  }

  QuirkMessenger quirkMessenger(quirkConfig());
  QSimMessenger qsimMessenger;
  if (seed > 0)
    G4Random::setTheSeed(seed);

  // the master never tracks: the quirk code relies on worker threads
  G4MTRunManager* rm = mt ? new G4MTRunManager : new G4TaskRunManager;
  rm->SetNumberOfThreads(nthreads);

  // phase 1: configuration lines
  std::vector<std::string> lines;
  {
    std::ifstream in(macro);
    if (!in) {
      G4cerr << "cannot open " << macro << G4endl;
      return 1;
    }
    for (std::string l; std::getline(in, l);)
      lines.push_back(l);
  }
  G4UImanager* ui = G4UImanager::GetUIpointer();
  size_t iline = 0;
  for (; iline < lines.size(); ++iline) {
    std::string l = lines[iline];
    l.erase(0, l.find_first_not_of(" \t"));
    if (l.empty() || l[0] == '#')
      continue;
    if (l.rfind("/run/initialize", 0) == 0)
      break;
    if (ui->ApplyCommand(l) != 0) {
      G4cerr << "command failed: " << l << G4endl;
      return 1;
    }
  }
  if (!output.empty())
    qsimConfig().output = output;

  const QSimConfig& c = qsimConfig();
  rm->SetUserInitialization(new DetectorConstruction(c.card));
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = factory.GetReferencePhysList(c.physList);
  if (nullptr == phys) {
    G4cerr << "unknown physics list " << c.physList << G4endl;
    return 1;
  }
  phys->SetDefaultCutValue(c.prodCut);
  if (quirkConfig().mass > 0.)
    phys->RegisterPhysics(new QuirkPhysics(quirkConfig()));
  rm->SetUserInitialization(phys);
  quirkMessenger.Lock();
  rm->SetUserInitialization(new ActionInitialization);
  G4cout << "quirksim: " << G4Version << ", " << c.physList << ", " << (mt ? "MT" : "tasking") << " run manager, "
         << rm->GetNumberOfThreads() << " threads" << G4endl;

  // phase 2: /run/initialize and the rest
  if (iline >= lines.size())
    ui->ApplyCommand("/run/initialize");
  for (; iline < lines.size(); ++iline) {
    std::string l = lines[iline];
    l.erase(0, l.find_first_not_of(" \t"));
    if (l.empty() || l[0] == '#')
      continue;
    if (ui->ApplyCommand(l) != 0) {
      G4cerr << "command failed: " << l << G4endl;
      return 1;
    }
  }
  if (nevents > 0)
    rm->BeamOn(nevents);
  delete rm;
  return 0;
}
