#ifndef GAGG_SIPM_DIGITIZER_HH
#define GAGG_SIPM_DIGITIZER_HH
#include "SiPMResponse.hh"
#include "G4VDigitizerModule.hh"
#include "G4VDigi.hh"
#include "G4Allocator.hh"
#include "G4TDigiCollection.hh"
struct WaveformEventData;

class SiPMDigi final : public G4VDigi {
public:
  sipm::Result response;
  void* operator new(std::size_t);
  void operator delete(void*);
};
using SiPMDigiCollection = G4TDigiCollection<SiPMDigi>;

class SiPMDigitizer final : public G4VDigitizerModule {
public:
  static constexpr const char* Name = "SiPMDigitizer";
  static constexpr const char* Collection = "SiPMDigitizer/SiPMWaveforms";
  SiPMDigitizer();
  void SetInput(const WaveformEventData& event) { input=&event; }
  void Digitize() override;
private:
  const WaveformEventData* input = nullptr;
};
#endif
