#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "Constants.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include <chrono>
#include <random>
#include <string>
#include <sys/stat.h>

#include "Randomize.hh"
#include "Rtypes.h"

namespace
{
  ULong64_t MakeRandomSeed()
  {
    const auto now = static_cast<ULong64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );

    std::random_device rd;
    const auto r0 = static_cast<ULong64_t>(rd());
    const auto r1 = static_cast<ULong64_t>(rd());

    // Mix clock entropy and random_device entropy.
    // The final seed is kept in the positive 31-bit range accepted by CLHEP/Geant4.
    ULong64_t mixed = now;
    mixed ^= (r0 << 1);
    mixed ^= (r1 << 33);
    mixed ^= (mixed >> 29);
    mixed *= 0x9E3779B97F4A7C15ULL;
    mixed ^= (mixed >> 32);

    return 1ULL + (mixed % 2147483646ULL);
  }

  void MakeDataDirectory()
  {
    struct stat st;
    if(stat(DATAPATH, &st) != 0){
      mkdir(DATAPATH, 0755);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc, char** argv)
{
  MakeDataDirectory();

  G4UIExecutive* ui = nullptr;
  if(argc == 1){
    ui = new G4UIExecutive(argc, argv);
  }

  auto run_manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  auto detector = new DetectorConstruction();
  run_manager->SetUserInitialization(detector);
  run_manager->SetUserInitialization(new PhysicsList());

  const auto random_seed = MakeRandomSeed();
  CLHEP::HepRandom::setTheSeed(static_cast<long>(random_seed));

  G4cout << "\n----> Random seed = " << random_seed << G4endl;

  G4String macro_file_name = "interactive";
  if(argc > 1){
    macro_file_name = argv[1];
  }

  run_manager->SetUserInitialization(new ActionInitialization(detector, random_seed, macro_file_name));

  auto vis_manager = new G4VisExecutive();
  vis_manager->Initialize();

  auto ui_manager = G4UImanager::GetUIpointer();

  if(!ui){
    G4String command = "/control/execute ";
    G4String file_name = argv[1];
    ui_manager->ApplyCommand(command + file_name);
  }else{
    ui_manager->ApplyCommand("/control/execute ../macros/init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete vis_manager;
  delete run_manager;

  return 0;
}
