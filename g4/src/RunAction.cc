#include "RunAction.hh"
#include "RootIO.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction(RootIO* io, DetectorConstruction* det)
: root_io(io), detector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
  if(root_io){
    delete root_io;
    root_io = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run*)
{
  if(root_io){
    if(detector){
      root_io->SetStepLimits(detector->GetGaggMaxStep(),
                             detector->GetMylarFrontMaxStep(),
                             detector->GetMylarSideMaxStep());
    }
    root_io->OpenFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run* run)
{
  if(root_io){
    if(run){
      root_io->SetNEvents(run->GetNumberOfEvent());
    }
    root_io->FillRunInfoTree();
    root_io->CloseFile();
  }
}
