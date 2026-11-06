#ifndef QSIM_QSimMessenger_hh
#define QSIM_QSimMessenger_hh

// /qsim/ UI commands filling QSimConfig (PreInit and Idle, not broadcast).

#include "G4UImessenger.hh"

#include <memory>
#include <vector>

class G4UIdirectory;
class G4UIcommand;

class QSimMessenger : public G4UImessenger {
public:
  QSimMessenger();
  ~QSimMessenger() override;
  void SetNewValue(G4UIcommand* cmd, G4String value) override;

private:
  std::vector<std::unique_ptr<G4UIdirectory>> m_dirs;
  std::vector<std::unique_ptr<G4UIcommand>> m_cmds;
};

#endif
