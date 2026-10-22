#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "RootIO.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SiPMDigitizer.hh"
#include "G4DigiManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ActionInitialization::ActionInitialization(DetectorConstruction* det, ULong64_t seed, const G4String& macro)
: detector(det), random_seed(seed), macro_file_name(macro)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::Build() const
{
  G4DigiManager::GetDMpointer()->AddNewModule(new SiPMDigitizer());
  auto root_io = new RootIO();
  root_io->SetRandomSeed(random_seed);
  root_io->SetMacroFile(macro_file_name);

  auto run_action = new RunAction(root_io, detector);
  SetUserAction(run_action);

  auto event_action = new EventAction(root_io);
  SetUserAction(event_action);

  SetUserAction(new PrimaryGeneratorAction(event_action));
  SetUserAction(new TrackingAction(event_action));
  SetUserAction(new SteppingAction(event_action, detector));
}
