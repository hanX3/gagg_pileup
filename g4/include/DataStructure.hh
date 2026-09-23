#ifndef DataStructure_H
#define DataStructure_H 1

#include <vector>
#include <string>
#include <cstdint>

#include "Constants.hh"

#include "Rtypes.h"
#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Internal particle category labels. ROOT output uses string fields instead of
// these numeric labels.
enum PrimaryType
{
  kPrimaryGamma   = 0,
  kPrimaryXray    = 1,
  kPrimaryAlpha   = 2,
  kPrimaryElectron= 3,
  kPrimaryOther   = 99
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Internal container used to keep one optical-photon SiPM hit together while
// sorting by arrival time. ROOT output is written to the separated vector
// branches: *_t_ps, *_x_mm, and *_y_mm.
struct OpticalPhotonHitData
{
  UInt_t t_ps;
  Float_t x_mm;
  Float_t y_mm;

  void Clear()
  {
    t_ps = 0;
    x_mm = 0.;
    y_mm = 0.;
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct RunInfoData
{
  UInt_t t_length_ps;
  UInt_t n_events;
  ULong64_t random_seed;
  Float_t gagg_max_step_um;
  Float_t mylar_front_max_step_um;
  Float_t mylar_side_max_step_um;
  TString macro_file;
  TString source_config_tag;

  void Clear()
  {
    t_length_ps = EVENT_TIME_LENGTH_PS;
    n_events = 0;
    random_seed = 0;
    gagg_max_step_um = DEFAULT_GAGG_MAX_STEP_UM;
    mylar_front_max_step_um = DEFAULT_MYLAR_FRONT_MAX_STEP_UM;
    mylar_side_max_step_um = DEFAULT_MYLAR_SIDE_MAX_STEP_UM;
    macro_file = "";
    source_config_tag = "";
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct WaveformEventData
{
  UInt_t event_id;
  UInt_t n_primary;
  ULong64_t n_optical_photons_arrived_total;
  Float_t edep_total_MeV;
  std::vector<OpticalPhotonHitData> all_photon_arrival_hits;
  std::vector<UInt_t> all_photon_arrival_t_ps;
  std::vector<Float_t> all_photon_arrival_x_mm;
  std::vector<Float_t> all_photon_arrival_y_mm;

  void Clear()
  {
    event_id = 0;
    n_primary = 0;
    n_optical_photons_arrived_total = 0;
    edep_total_MeV = 0.;
    all_photon_arrival_hits.clear();
    all_photon_arrival_t_ps.clear();
    all_photon_arrival_x_mm.clear();
    all_photon_arrival_y_mm.clear();
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct PrimaryPhotonData
{
  UInt_t event_id;
  UInt_t primary_id;
  std::string primary_particle_name;
  std::string primary_response_name;
  Float_t primary_energy_MeV;
  UInt_t primary_time_ps;
  Float_t primary_x_mm;
  Float_t primary_y_mm;
  Float_t primary_z_mm;
  Float_t primary_px;
  Float_t primary_py;
  Float_t primary_pz;

  Float_t edep_total_MeV;
  UInt_t n_scint_photons_generated;
  UInt_t n_optical_photons_arrived;
  std::vector<OpticalPhotonHitData> photon_arrival_hits;
  std::vector<UInt_t> photon_arrival_t_ps;
  std::vector<Float_t> photon_arrival_x_mm;
  std::vector<Float_t> photon_arrival_y_mm;

  void Clear()
  {
    event_id = 0;
    primary_id = 0;
    primary_particle_name.clear();
    primary_response_name.clear();
    primary_energy_MeV = 0.;
    primary_time_ps = 0;
    primary_x_mm = 0.;
    primary_y_mm = 0.;
    primary_z_mm = 0.;
    primary_px = 0.;
    primary_py = 0.;
    primary_pz = 1.;

    edep_total_MeV = 0.;
    n_scint_photons_generated = 0;
    n_optical_photons_arrived = 0;
    photon_arrival_hits.clear();
    photon_arrival_t_ps.clear();
    photon_arrival_x_mm.clear();
    photon_arrival_y_mm.clear();
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
