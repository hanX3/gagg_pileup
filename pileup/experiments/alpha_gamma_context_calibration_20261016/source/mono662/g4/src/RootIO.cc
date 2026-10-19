#include "RootIO.hh"
#include "RunConfig.hh"
#include "SiPMConfig.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <thread>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

namespace
{
  std::string TrimLeft(const std::string& s)
  {
    const auto it = std::find_if(s.begin(), s.end(), [](unsigned char ch){
      return !std::isspace(ch);
    });
    return std::string(it, s.end());
  }

  bool IsSourceEnableOrRateLine(const std::string& line)
  {
    const auto trimmed = TrimLeft(line);
    if(trimmed.empty()) return false;
    if(trimmed[0] == '#') return false;

    const bool is_source_line = (trimmed.find("/gagg/source/") == 0);
    const bool is_enable_or_rate =
      (trimmed.find("/enable") != std::string::npos) ||
      (trimmed.find("/rateHz") != std::string::npos);

    return is_source_line && is_enable_or_rate;
  }

  std::string NormalizeSourceConfigLine(const std::string& line)
  {
    std::string s = TrimLeft(line);

    const std::string prefix = "/gagg/source/";
    if(s.find(prefix) == 0){
      s.erase(0, prefix.size());
    }

    std::string out;
    out.reserve(s.size());

    bool last_was_separator = false;
    for(const auto ch_raw : s){
      const unsigned char ch = static_cast<unsigned char>(ch_raw);
      const bool is_separator = (ch_raw == '/') || std::isspace(ch);

      if(is_separator){
        if(!out.empty() && !last_was_separator){
          out.push_back('_');
          last_was_separator = true;
        }
      }else{
        out.push_back(ch_raw);
        last_was_separator = false;
      }
    }

    while(!out.empty() && out.back() == '_'){
      out.pop_back();
    }

    return out;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO::RootIO()
{
  run_info_data.Clear();
  waveform_event_data.Clear();
  primary_photon_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO::~RootIO()
{
  if(root_file){
    delete root_file;
    root_file = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::SetNEvents(G4int n_events)
{
  run_info_data.n_events = static_cast<UInt_t>(n_events);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::SetRandomSeed(ULong64_t seed)
{
  run_info_data.random_seed = seed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::SetStepLimits(G4double gagg_step, G4double mylar_front_step, G4double mylar_side_step)
{
  run_info_data.gagg_max_step_um = static_cast<Float_t>(gagg_step/um);
  run_info_data.mylar_front_max_step_um = static_cast<Float_t>(mylar_front_step/um);
  run_info_data.mylar_side_max_step_um = static_cast<Float_t>(mylar_side_step/um);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::SetMacroFile(const G4String& macro_file_name)
{
  this->macro_file_name = macro_file_name;
  run_info_data.macro_file = macro_file_name.c_str();

  if(macro_file_name.empty() || macro_file_name == "interactive"){
    run_info_data.source_config_tag = "";
    return;
  }

  std::ifstream fin(macro_file_name.c_str());
  if(!fin){
    run_info_data.source_config_tag = "ERROR_cannot_open_macro_file";
    return;
  }

  std::ostringstream tag;
  std::string line;
  while(std::getline(fin, line)){
    run_info_data.macro_contents += (line + "\n").c_str();
    if(IsSourceEnableOrRateLine(line)){
      const auto item = NormalizeSourceConfigLine(line);
      if(item.empty()) continue;

      if(tag.tellp() > 0){
        tag << "_";
      }
      tag << item;
    }
  }

  run_info_data.source_config_tag = tag.str().c_str();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::OpenFile()
{
  const auto& config = RunConfig::Instance();
  config.Validate();
  auto& sensor = SiPMConfig::Instance();
  sensor.Prepare(config.WindowPS()*0.001);
  sipm_enabled = sensor.Enabled();
  sipm_config_file = sensor.ResolvedFile();
  sipm_config_json = sipm_enabled ? sensor.Model().json_text : "";
  sipm_avalanche_base_seed = sipm_enabled ? sensor.AvalancheSeed() : 0;
  sipm_readout_base_seed = sipm_enabled ? sensor.ReadoutSeed() : 0;
  sipm_model_version = "gagg-native-sipm-v1";
  sipm_rng = "std::mt19937_64; splitmix64 event/substream seeds; C++ standard distributions; independent of CLHEP";
  run_info_data.t_length_ps = config.WindowPS();
  run_info_data.analysis_start_ps = config.PreWindowPS();
  run_info_data.analysis_end_ps = config.WindowPS() - config.PostWindowPS();
  run_info_data.random_seed = static_cast<ULong64_t>(CLHEP::HepRandom::getTheSeed());
  std::ostringstream engine_state;
  CLHEP::HepRandom::saveFullState(engine_state);
  run_info_data.random_engine_state = engine_state.str().c_str();
  const auto* source = dynamic_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if(source){ source->FillRunInfo(run_info_data); }

  std::filesystem::path path;
  if(config.OutputFile().empty()){
    do {
      const auto t = time(nullptr);
      const auto* tt = localtime(&t);
      strftime(file_name, sizeof(file_name), "%Y%m%d_%Hh%Mm%Ss", tt);
      path = std::filesystem::absolute(std::filesystem::path(DATAPATH) /
        (std::string(file_name) + ".root"));
      if(!std::filesystem::exists(path)) break;
      std::this_thread::sleep_for(std::chrono::seconds(1));
    } while(true);
  }else{
    path = std::filesystem::absolute(config.OutputFile().c_str());
  }
  std::filesystem::create_directories(path.parent_path());
  if(std::filesystem::exists(path)){
    throw std::runtime_error("ROOT output already exists: " + path.string());
  }
  if(root_file){ delete root_file; root_file = nullptr; }
  root_file = TFile::Open(path.string().c_str(), "NEW");
  if(!root_file || root_file->IsZombie()){
    throw std::runtime_error("Cannot create ROOT output: " + path.string());
  }
  G4cout << "----> ROOT output: " << path.string() << G4endl;

  // run info
  run_info_tree = new TTree("RunInfo", "run information");
  run_info_tree->Branch("t_length_ps", &run_info_data.t_length_ps, "t_length_ps/i");
  run_info_tree->Branch("n_events", &run_info_data.n_events, "n_events/i");
  run_info_tree->Branch("random_seed", &run_info_data.random_seed, "random_seed/l");
  run_info_tree->Branch("gagg_max_step_um", &run_info_data.gagg_max_step_um, "gagg_max_step_um/F");
  run_info_tree->Branch("mylar_front_max_step_um", &run_info_data.mylar_front_max_step_um, "mylar_front_max_step_um/F");
  run_info_tree->Branch("mylar_side_max_step_um", &run_info_data.mylar_side_max_step_um, "mylar_side_max_step_um/F");
  run_info_tree->Branch("macro_file", &run_info_data.macro_file);
  run_info_tree->Branch("source_config_tag", &run_info_data.source_config_tag);
  run_info_tree->Branch("macro_contents", &run_info_data.macro_contents);
  run_info_tree->Branch("random_engine_state", &run_info_data.random_engine_state);
  run_info_tree->Branch("analysis_start_ps", &run_info_data.analysis_start_ps, "analysis_start_ps/i");
  run_info_tree->Branch("analysis_end_ps", &run_info_data.analysis_end_ps, "analysis_end_ps/i");
  run_info_tree->Branch("pileup_enabled", &run_info_data.pileup_enabled, "pileup_enabled/O");
  run_info_tree->Branch("brems_enabled", &run_info_data.brems_enabled, "brems_enabled/O");
  run_info_tree->Branch("c12_capture_enabled", &run_info_data.c12_capture_enabled, "c12_capture_enabled/O");
  run_info_tree->Branch("alpha_p11b_enabled", &run_info_data.alpha_p11b_enabled, "alpha_p11b_enabled/O");
  run_info_tree->Branch("brems_rate_hz", &run_info_data.brems_rate_hz, "brems_rate_hz/D");
  run_info_tree->Branch("c12_capture_rate_hz", &run_info_data.c12_capture_rate_hz, "c12_capture_rate_hz/D");
  run_info_tree->Branch("alpha_p11b_rate_hz", &run_info_data.alpha_p11b_rate_hz, "alpha_p11b_rate_hz/D");
  run_info_tree->Branch("use_disk_cone_source", &run_info_data.use_disk_cone_source, "use_disk_cone_source/O");
  run_info_tree->Branch("aim_at_gagg", &run_info_data.aim_at_gagg, "aim_at_gagg/O");
  run_info_tree->Branch("source_disk_radius_mm", &run_info_data.source_disk_radius_mm, "source_disk_radius_mm/D");
  run_info_tree->Branch("source_cone_half_angle_deg", &run_info_data.source_cone_half_angle_deg, "source_cone_half_angle_deg/D");
  run_info_tree->Branch("sipm_enabled", &sipm_enabled, "sipm_enabled/O");
  run_info_tree->Branch("sipm_schema_version", &sipm_schema_version, "sipm_schema_version/i");
  run_info_tree->Branch("sipm_avalanche_base_seed", &sipm_avalanche_base_seed, "sipm_avalanche_base_seed/l");
  run_info_tree->Branch("sipm_readout_base_seed", &sipm_readout_base_seed, "sipm_readout_base_seed/l");
  run_info_tree->Branch("sipm_config_file", &sipm_config_file);
  run_info_tree->Branch("sipm_config_json", &sipm_config_json);
  run_info_tree->Branch("sipm_model_version", &sipm_model_version);
  run_info_tree->Branch("sipm_rng", &sipm_rng);

  sipm_tree = nullptr;
  if(sipm_enabled){
    sipm_tree = new TTree("SiPMEvent", "One G4 digitized waveform per event; device parameters and calibration status in RunInfo");
#define SIPM_SCALAR(field, type) sipm_tree->Branch(#field, &sipm_data.field, #field "/" type)
    SIPM_SCALAR(event_id,"i"); SIPM_SCALAR(avalanche_seed,"l"); SIPM_SCALAR(readout_seed,"l");
    SIPM_SCALAR(n_input_photons,"l"); SIPM_SCALAR(n_outside_sensor,"l");
    SIPM_SCALAR(n_delayed_photons_outside,"l"); SIPM_SCALAR(n_secondary_outside,"l"); SIPM_SCALAR(n_candidates,"l");
    SIPM_SCALAR(n_photon_avalanches,"l"); SIPM_SCALAR(n_dark_avalanches,"l");
    SIPM_SCALAR(n_crosstalk_avalanches,"l"); SIPM_SCALAR(n_afterpulse_avalanches,"l"); SIPM_SCALAR(n_clipped_samples,"l");
    SIPM_SCALAR(charge_total_C,"D"); SIPM_SCALAR(charge_gate_C,"D"); SIPM_SCALAR(current_integral_C,"D");
    SIPM_SCALAR(current_tail_C,"D"); SIPM_SCALAR(charge_balance_error_C,"D");
#undef SIPM_SCALAR
#define SIPM_VECTOR(field) sipm_tree->Branch(#field, &sipm_data.field)
    SIPM_VECTOR(avalanche_time_ns); SIPM_VECTOR(avalanche_charge_C); SIPM_VECTOR(avalanche_recovery);
    SIPM_VECTOR(avalanche_cell); SIPM_VECTOR(avalanche_kind); SIPM_VECTOR(avalanche_parent); SIPM_VECTOR(avalanche_photon_index);
    SIPM_VECTOR(current_time_ns); SIPM_VECTOR(current_A); SIPM_VECTOR(sample_time_ns);
    SIPM_VECTOR(signal_mV); SIPM_VECTOR(analog_mV); SIPM_VECTOR(adc_counts); SIPM_VECTOR(adc_clipped);
#undef SIPM_VECTOR
  }

  // waveform event
  waveform_event_tree = new TTree("WaveformEvent", "one entry is one independent waveform window; duration in RunInfo");
  waveform_event_tree->Branch("event_id", &waveform_event_data.event_id, "event_id/i");
  waveform_event_tree->Branch("n_primary", &waveform_event_data.n_primary, "n_primary/i");
  waveform_event_tree->Branch("n_optical_photons_arrived_total", &waveform_event_data.n_optical_photons_arrived_total, "n_optical_photons_arrived_total/l");
  waveform_event_tree->Branch("edep_total_MeV", &waveform_event_data.edep_total_MeV, "edep_total_MeV/F");
  waveform_event_tree->Branch("all_photon_arrival_t_ps", &waveform_event_data.all_photon_arrival_t_ps);
  waveform_event_tree->Branch("all_photon_arrival_x_mm", &waveform_event_data.all_photon_arrival_x_mm);
  waveform_event_tree->Branch("all_photon_arrival_y_mm", &waveform_event_data.all_photon_arrival_y_mm);
  waveform_event_tree->Branch("photon_energy_eV", &waveform_event_data.photon_energy_eV);

  // primary photon
  primary_photon_tree = new TTree("PrimaryPhoton", "one entry is one primary response");
  primary_photon_tree->Branch("event_id", &primary_photon_data.event_id, "event_id/i");
  primary_photon_tree->Branch("primary_id", &primary_photon_data.primary_id, "primary_id/i");
  primary_photon_tree->Branch("primary_particle_name", &primary_photon_data.primary_particle_name);
  primary_photon_tree->Branch("primary_response_name", &primary_photon_data.primary_response_name);
  primary_photon_tree->Branch("primary_energy_MeV", &primary_photon_data.primary_energy_MeV, "primary_energy_MeV/F");
  primary_photon_tree->Branch("primary_time_ps", &primary_photon_data.primary_time_ps, "primary_time_ps/i");
  primary_photon_tree->Branch("primary_x_mm", &primary_photon_data.primary_x_mm, "primary_x_mm/F");
  primary_photon_tree->Branch("primary_y_mm", &primary_photon_data.primary_y_mm, "primary_y_mm/F");
  primary_photon_tree->Branch("primary_z_mm", &primary_photon_data.primary_z_mm, "primary_z_mm/F");
  primary_photon_tree->Branch("primary_px", &primary_photon_data.primary_px, "primary_px/F");
  primary_photon_tree->Branch("primary_py", &primary_photon_data.primary_py, "primary_py/F");
  primary_photon_tree->Branch("primary_pz", &primary_photon_data.primary_pz, "primary_pz/F");
  primary_photon_tree->Branch("edep_total_MeV", &primary_photon_data.edep_total_MeV, "edep_total_MeV/F");
  primary_photon_tree->Branch("n_scint_photons_generated", &primary_photon_data.n_scint_photons_generated, "n_scint_photons_generated/i");
  primary_photon_tree->Branch("n_optical_photons_arrived", &primary_photon_data.n_optical_photons_arrived, "n_optical_photons_arrived/i");
  primary_photon_tree->Branch("photon_arrival_t_ps", &primary_photon_data.photon_arrival_t_ps);
  primary_photon_tree->Branch("photon_arrival_x_mm", &primary_photon_data.photon_arrival_x_mm);
  primary_photon_tree->Branch("photon_arrival_y_mm", &primary_photon_data.photon_arrival_y_mm);
  primary_photon_tree->Branch("photon_energy_eV", &primary_photon_data.photon_energy_eV);

  if(!run_info_tree || !waveform_event_tree || !primary_photon_tree){
    G4cout << "\n can't create ROOT trees" << G4endl;
    return;
  }

  G4cout << "\n----> ROOT trees are created in " << path.string() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillRunInfoTree()
{
  if(run_info_tree){
    run_info_tree->Fill();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillWaveformEventTree(WaveformEventData& data)
{
  waveform_event_data = data;
  waveform_event_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillPrimaryPhotonTree(PrimaryPhotonData& data)
{
  primary_photon_data = data;
  primary_photon_tree->Fill();
}

void RootIO::FillSiPMEventTree(const sipm::Result& data)
{
  if(!sipm_tree) throw std::runtime_error("SiPMEvent tree was not prepared for this run.");
  sipm_data = data;
  sipm_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::CloseFile()
{
  if(!root_file){
    return;
  }

  root_file->cd();
  if(run_info_tree){ run_info_tree->Write(); }
  if(waveform_event_tree){ waveform_event_tree->Write(); }
  if(primary_photon_tree){ primary_photon_tree->Write(); }
  if(sipm_tree){ sipm_tree->Write(); }
  root_file->Close();

  G4cout << "\n----> ROOT file is saved.\n" << G4endl;
}
