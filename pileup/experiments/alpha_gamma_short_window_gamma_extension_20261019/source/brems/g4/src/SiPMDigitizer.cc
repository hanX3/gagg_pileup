#include "SiPMDigitizer.hh"
#include "SiPMConfig.hh"
#include "DataStructure.hh"
#include "RunConfig.hh"
#include <memory>
#include <stdexcept>

namespace { G4ThreadLocal G4Allocator<SiPMDigi>* allocator=nullptr; }
void* SiPMDigi::operator new(std::size_t) {
  if(!allocator) allocator=new G4Allocator<SiPMDigi>;
  return allocator->MallocSingle();
}
void SiPMDigi::operator delete(void* object) { allocator->FreeSingle(static_cast<SiPMDigi*>(object)); }

SiPMDigitizer::SiPMDigitizer() : G4VDigitizerModule(Name) { collectionName.push_back("SiPMWaveforms"); }
void SiPMDigitizer::Digitize() {
  if(!input) throw std::runtime_error("SiPMDigitizer has no current waveform input.");
  const auto* event=input; input=nullptr;
  const auto& settings=SiPMConfig::Instance();
  const auto& run=RunConfig::Instance();
  std::vector<sipm::Photon> photons;
  photons.reserve(event->all_photon_arrival_hits.size());
  for(const auto& p:event->all_photon_arrival_hits)
    photons.push_back({p.t_ps*0.001,p.x_mm,p.y_mm,p.photon_energy_eV});
  auto digi=std::make_unique<SiPMDigi>();
  digi->response=sipm::Simulate(settings.Model(),photons,run.WindowPS()*0.001,
    run.PreWindowPS()*0.001,(run.WindowPS()-run.PostWindowPS())*0.001,
    sipm::EventSeed(settings.AvalancheSeed(),event->event_id),sipm::EventSeed(settings.ReadoutSeed(),event->event_id));
  digi->response.event_id=event->event_id;
  auto collection=std::make_unique<SiPMDigiCollection>(moduleName,collectionName[0]);
  collection->insert(digi.release());
  StoreDigiCollection(collection.release()); // Ownership transfers to this G4Event.
}
