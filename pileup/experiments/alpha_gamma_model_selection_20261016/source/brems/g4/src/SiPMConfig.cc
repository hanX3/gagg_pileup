#include "SiPMConfig.hh"
#include "G4GenericMessenger.hh"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

SiPMConfig& SiPMConfig::Instance() {
  static auto* settings = new SiPMConfig;
  return *settings;
}
SiPMConfig::SiPMConfig() {
  messenger=new G4GenericMessenger(this,"/gagg/sipm/","Native SiPM digitization; all configuration units are explicit in JSON.");
  messenger->DeclareProperty("enable",enabled).SetStates(G4State_PreInit,G4State_Idle);
  messenger->DeclareProperty("configFile",config_file,"JSON path relative to the working directory.")
    .SetStates(G4State_PreInit,G4State_Idle);
  messenger->DeclareProperty("seed",avalanche_seed,"Independent avalanche RNG base seed.")
    .SetParameterName("seed",false).SetRange("seed>0 && seed<2147483647").SetStates(G4State_PreInit,G4State_Idle);
  messenger->DeclareProperty("readoutSeed",readout_seed,"Independent electronics RNG base seed.")
    .SetParameterName("seed",false).SetRange("seed>0 && seed<2147483647").SetStates(G4State_PreInit,G4State_Idle);
}
void SiPMConfig::Prepare(double window_ns) {
  model.reset(); resolved_file.clear();
  if(!enabled) return;
  if(config_file.empty()) throw std::runtime_error("SiPM enabled: /gagg/sipm/configFile is required.");
  resolved_file=std::filesystem::absolute(config_file.c_str()).lexically_normal().string();
  std::ifstream input(resolved_file);
  if(!input) throw std::runtime_error("Cannot open SiPM config: "+resolved_file);
  std::ostringstream contents; contents<<input.rdbuf();
  auto candidate=std::make_unique<sipm::Config>(sipm::Config::FromJSON(contents.str()));
  candidate->ValidateWindow(window_ns);
  model=std::move(candidate);
}
const sipm::Config& SiPMConfig::Model() const {
  if(!model) throw std::runtime_error("SiPM configuration has not been prepared for this run.");
  return *model;
}
