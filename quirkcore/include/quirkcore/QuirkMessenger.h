#ifndef QUIRKCORE_QuirkMessenger_H
#define QUIRKCORE_QuirkMessenger_H

// /quirk/ UI commands filling a QuirkConfig. The particle definitions are made
// when the physics list is registered, so the commands only act in PreInit,
// before that (quirksim applies the /quirk/ lines of the macro first).

#include "quirkcore/QuirkConfig.h"

#include "G4UImessenger.hh"

#include <memory>
#include <vector>

class G4UIdirectory;
class G4UIcommand;

class QuirkMessenger : public G4UImessenger {
public:
  explicit QuirkMessenger(QuirkConfig& cfg);
  ~QuirkMessenger() override;

  void SetNewValue(G4UIcommand* cmd, G4String value) override;
  G4String GetCurrentValue(G4UIcommand* cmd) override;

  // after this the particles exist: further changes are refused
  void Lock() { m_locked = true; }

private:
  QuirkConfig& m_cfg;
  G4bool m_locked{false};
  std::unique_ptr<G4UIdirectory> m_dir;
  std::vector<std::unique_ptr<G4UIcommand>> m_cmds;
};

#endif
