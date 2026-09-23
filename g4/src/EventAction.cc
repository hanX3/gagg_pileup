#include "EventAction.hh"
#include "RootIO.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

namespace
{
  void SortOpticalPhotonHits(std::vector<OpticalPhotonHitData>& hits)
  {
    std::stable_sort(hits.begin(), hits.end(),
      [](const OpticalPhotonHitData& a, const OpticalPhotonHitData& b){
        return a.t_ps < b.t_ps;
      });
  }

  void FillOpticalPhotonHitBranches(const std::vector<OpticalPhotonHitData>& hits,
                                    std::vector<UInt_t>& t_ps,
                                    std::vector<Float_t>& x_mm,
                                    std::vector<Float_t>& y_mm)
  {
    t_ps.clear();
    x_mm.clear();
    y_mm.clear();

    t_ps.reserve(hits.size());
    x_mm.reserve(hits.size());
    y_mm.reserve(hits.size());

    for(const auto& hit : hits){
      t_ps.push_back(hit.t_ps);
      x_mm.push_back(hit.x_mm);
      y_mm.push_back(hit.y_mm);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction(RootIO* io)
: root_io(io)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event* event)
{
  waveform_event_data.Clear();
  primary_photon_data_vec.clear();

  waveform_event_data.event_id = event->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event*)
{
  waveform_event_data.n_primary = primary_photon_data_vec.size();

  SortOpticalPhotonHits(waveform_event_data.all_photon_arrival_hits);
  FillOpticalPhotonHitBranches(waveform_event_data.all_photon_arrival_hits,
                               waveform_event_data.all_photon_arrival_t_ps,
                               waveform_event_data.all_photon_arrival_x_mm,
                               waveform_event_data.all_photon_arrival_y_mm);
  waveform_event_data.n_optical_photons_arrived_total = waveform_event_data.all_photon_arrival_t_ps.size();

  if(root_io){
    root_io->FillWaveformEventTree(waveform_event_data);
  }

  for(auto& primary_data : primary_photon_data_vec){
    SortOpticalPhotonHits(primary_data.photon_arrival_hits);
    FillOpticalPhotonHitBranches(primary_data.photon_arrival_hits,
                                 primary_data.photon_arrival_t_ps,
                                 primary_data.photon_arrival_x_mm,
                                 primary_data.photon_arrival_y_mm);
    primary_data.n_optical_photons_arrived = primary_data.photon_arrival_t_ps.size();

    if(root_io){
      root_io->FillPrimaryPhotonTree(primary_data);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddPrimary(G4int primary_id,
                             const G4String& particle_name,
                             const G4String& response_name,
                             G4double energy,
                             G4double time,
                             const G4ThreeVector& position,
                             const G4ThreeVector& direction)
{
  PrimaryPhotonData data;
  data.Clear();

  data.event_id = waveform_event_data.event_id;
  data.primary_id = primary_id;
  data.primary_particle_name = particle_name;
  data.primary_response_name = response_name;
  data.primary_energy_MeV = energy / MeV;
  data.primary_time_ps = static_cast<UInt_t>(std::llround(time / ps));
  data.primary_x_mm = position.x() / mm;
  data.primary_y_mm = position.y() / mm;
  data.primary_z_mm = position.z() / mm;
  data.primary_px = direction.x();
  data.primary_py = direction.y();
  data.primary_pz = direction.z();

  if(primary_id >= static_cast<G4int>(primary_photon_data_vec.size())){
    primary_photon_data_vec.resize(primary_id + 1);
  }
  primary_photon_data_vec[primary_id] = data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool EventAction::IsValidPrimaryId(G4int primary_id) const
{
  if(primary_id < 0){
    return false;
  }
  if(primary_id >= static_cast<G4int>(primary_photon_data_vec.size())){
    return false;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddEdep(G4int primary_id, G4double edep)
{
  if(!IsValidPrimaryId(primary_id)){
    return;
  }

  G4double edep_MeV = edep / MeV;

  primary_photon_data_vec[primary_id].edep_total_MeV += edep_MeV;
  waveform_event_data.edep_total_MeV += edep_MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddScintPhotonGenerated(G4int primary_id)
{
  if(!IsValidPrimaryId(primary_id)){
    return;
  }

  primary_photon_data_vec[primary_id].n_scint_photons_generated++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddPhotonArrival(G4int primary_id, G4double time, const G4ThreeVector& position)
{
  if(!IsValidPrimaryId(primary_id)){
    return;
  }

  auto time_ps_64 = static_cast<unsigned long long>(std::llround(time / ps));
  if(time_ps_64 > EVENT_TIME_LENGTH_PS){
    return;
  }

  OpticalPhotonHitData hit;
  hit.t_ps = static_cast<UInt_t>(time_ps_64);
  hit.x_mm = static_cast<Float_t>(position.x() / mm);
  hit.y_mm = static_cast<Float_t>(position.y() / mm);

  waveform_event_data.all_photon_arrival_hits.push_back(hit);
  primary_photon_data_vec[primary_id].photon_arrival_hits.push_back(hit);
  primary_photon_data_vec[primary_id].n_optical_photons_arrived++;
}
