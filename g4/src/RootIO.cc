#include "RootIO.hh"

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cstring>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

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
void RootIO::OpenFile()
{
  time_t t;
  struct tm* tt;
  t=time(0);
  tt=localtime(&t);
  sprintf(file_name, "%d%02d%02d_%02dh%02dm%02ds", tt->tm_year+1900, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  G4cout << "\n----> ROOT file is opened in " << file_name << G4endl;

  std::stringstream ss;
  ss.str("");
  ss << DATAPATH << "/gagg_waveform_" << file_name << ".root";

  root_file = new TFile(ss.str().c_str(), "RECREATE");
  if(!root_file){
    G4cout << " RootIO::" << " problem creating the ROOT TFile!!!" << G4endl;
    return;
  }else{
    G4cout << " RootIO::" << " successful creating the " << ss.str().c_str() << "  !!!" << G4endl;
  }

  // run info
  run_info_tree = new TTree("RunInfo", "run information");
  run_info_tree->Branch("t_length_ps", &run_info_data.t_length_ps, "t_length_ps/i");
  run_info_tree->Branch("n_events", &run_info_data.n_events, "n_events/i");
  run_info_tree->Branch("random_seed", &run_info_data.random_seed, "random_seed/l");
  run_info_tree->Branch("gagg_max_step_um", &run_info_data.gagg_max_step_um, "gagg_max_step_um/F");
  run_info_tree->Branch("mylar_front_max_step_um", &run_info_data.mylar_front_max_step_um, "mylar_front_max_step_um/F");
  run_info_tree->Branch("mylar_side_max_step_um", &run_info_data.mylar_side_max_step_um, "mylar_side_max_step_um/F");

  // waveform event
  waveform_event_tree = new TTree("WaveformEvent", "one entry is one 1 ms waveform event");
  waveform_event_tree->Branch("event_id", &waveform_event_data.event_id, "event_id/i");
  waveform_event_tree->Branch("n_primary", &waveform_event_data.n_primary, "n_primary/i");
  waveform_event_tree->Branch("n_optical_photons_arrived_total", &waveform_event_data.n_optical_photons_arrived_total, "n_optical_photons_arrived_total/l");
  waveform_event_tree->Branch("edep_total_MeV", &waveform_event_data.edep_total_MeV, "edep_total_MeV/F");
  waveform_event_tree->Branch("all_photon_arrival_t_ps", &waveform_event_data.all_photon_arrival_t_ps);

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

  if(!run_info_tree || !waveform_event_tree || !primary_photon_tree){
    G4cout << "\n can't create ROOT trees" << G4endl;
    return;
  }

  G4cout << "\n----> ROOT trees are created in " << ss.str() << G4endl;
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
  root_file->Close();

  G4cout << "\n----> ROOT file is saved.\n" << G4endl;
}
