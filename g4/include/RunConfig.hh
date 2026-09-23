#ifndef GAGG_RUN_CONFIG_HH
#define GAGG_RUN_CONFIG_HH

#include "globals.hh"
#include <cstdint>

class G4GenericMessenger;

// Runtime settings shared by the serial source, recorder and RunInfo writer.
class RunConfig
{
public:
  static RunConfig& Instance();
  void Validate() const;
  std::uint32_t WindowPS() const;
  std::uint32_t PreWindowPS() const;
  std::uint32_t PostWindowPS() const;
  const G4String& OutputFile() const { return output_file; }

private:
  RunConfig();
  void SetSeed(G4int seed);
  G4GenericMessenger* messenger = nullptr;
  G4double window_length;
  G4double pre_window = 0.;
  G4double post_window = 0.;
  G4String output_file;
};

#endif
