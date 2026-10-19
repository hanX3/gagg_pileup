#ifndef RootIO_H
#define RootIO_H 1

#include "Constants.hh"
#include "DataStructure.hh"
#include "SiPMResponse.hh"
#include <globals.hh>

#include <fstream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RootIO
{
public:
  RootIO();
  ~RootIO();

// file
public:
  void OpenFile();
  void CloseFile();
  void SetNEvents(G4int n_events);
  void SetRandomSeed(ULong64_t seed);
  void SetStepLimits(G4double gagg_step, G4double mylar_front_step, G4double mylar_side_step);
  void SetMacroFile(const G4String& macro_file_name);

// run info
public:
  void FillRunInfoTree();

// waveform event
public:
  void FillWaveformEventTree(WaveformEventData& data);

// primary photon
public:
  void FillPrimaryPhotonTree(PrimaryPhotonData& data);
  void FillSiPMEventTree(const sipm::Result& data);

// private
private:

private:
  RunInfoData run_info_data;
  WaveformEventData waveform_event_data;
  PrimaryPhotonData primary_photon_data;
  sipm::Result sipm_data;
  Bool_t sipm_enabled = false;
  UInt_t sipm_schema_version = 1;
  ULong64_t sipm_avalanche_base_seed = 0, sipm_readout_base_seed = 0;
  std::string sipm_config_file, sipm_config_json, sipm_model_version, sipm_rng;

  G4String macro_file_name;

  TFile* root_file = nullptr;
  TTree* run_info_tree = nullptr;
  TTree* waveform_event_tree = nullptr;
  TTree* primary_photon_tree = nullptr;
  TTree* sipm_tree = nullptr;

  char file_name[1024];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
