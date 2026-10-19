#ifndef GAGG_SIPM_CONFIG_HH
#define GAGG_SIPM_CONFIG_HH
#include "SiPMResponse.hh"
#include "globals.hh"
#include <memory>
class G4GenericMessenger;
class SiPMConfig {
public:
  static SiPMConfig& Instance();
  bool Enabled() const { return enabled; }
  void Prepare(double window_ns);
  const sipm::Config& Model() const;
  const std::string& ResolvedFile() const { return resolved_file; }
  std::uint64_t AvalancheSeed() const { return avalanche_seed; }
  std::uint64_t ReadoutSeed() const { return readout_seed; }
private:
  SiPMConfig();
  G4GenericMessenger* messenger = nullptr;
  G4bool enabled = false;
  G4String config_file;
  G4int avalanche_seed = 1, readout_seed = 2;
  std::string resolved_file;
  std::unique_ptr<sipm::Config> model;
};
#endif
