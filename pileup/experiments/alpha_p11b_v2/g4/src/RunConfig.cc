#include "RunConfig.hh"
#include "Constants.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <limits>
#include <stdexcept>

RunConfig& RunConfig::Instance()
{
  // Process lifetime: UI commands must outlive all actions in the serial run.
  static auto* config = new RunConfig;
  return *config;
}

RunConfig::RunConfig() : window_length(EVENT_TIME_LENGTH_PS*ps)
{
  messenger = new G4GenericMessenger(this, "/gagg/run/", "Reproducible waveform runs.");
  auto& window = messenger->DeclarePropertyWithUnit(
    "windowLength", "us", window_length, "Full source and photon-recording interval [0,T).");
  window.SetParameterName("length", false).SetRange("length>0")
    .SetStates(G4State_PreInit, G4State_Idle);
  auto& pre = messenger->DeclarePropertyWithUnit(
    "preWindow", "us", pre_window, "History interval before the central analysis gate.");
  pre.SetParameterName("length", false).SetRange("length>=0")
    .SetStates(G4State_PreInit, G4State_Idle);
  auto& post = messenger->DeclarePropertyWithUnit(
    "postWindow", "us", post_window, "Response guard interval after the central analysis gate.");
  post.SetParameterName("length", false).SetRange("length>=0")
    .SetStates(G4State_PreInit, G4State_Idle);
  messenger->DeclareProperty("outputFile", output_file,
    "ROOT path relative to the working directory; an existing file is never overwritten.")
    .SetStates(G4State_PreInit, G4State_Idle);
  messenger->DeclareMethod("seed", &RunConfig::SetSeed, "Set a reproducible CLHEP seed.")
    .SetParameterName("seed", false).SetRange("seed>0 && seed<2147483647")
    .SetStates(G4State_PreInit, G4State_Idle);
}

void RunConfig::SetSeed(G4int seed)
{
  CLHEP::HepRandom::setTheSeed(seed);
}

void RunConfig::Validate() const
{
  const auto max_ps = static_cast<double>(std::numeric_limits<std::uint32_t>::max());
  for(const auto value : {window_length, pre_window, post_window}){
    if(!std::isfinite(value) || value < 0. || value/ps > max_ps){
      throw std::runtime_error("Run window must be finite, nonnegative and fit uint32 picoseconds.");
    }
  }
  if(WindowPS() == 0 || static_cast<std::uint64_t>(PreWindowPS()) + PostWindowPS() >= WindowPS()){
    throw std::runtime_error("windowLength must exceed preWindow + postWindow at picosecond precision.");
  }
}

std::uint32_t RunConfig::WindowPS() const
{ return static_cast<std::uint32_t>(std::llround(window_length/ps)); }
std::uint32_t RunConfig::PreWindowPS() const
{ return static_cast<std::uint32_t>(std::llround(pre_window/ps)); }
std::uint32_t RunConfig::PostWindowPS() const
{ return static_cast<std::uint32_t>(std::llround(post_window/ps)); }
