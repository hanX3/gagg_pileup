#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "TrackInformation.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* action, DetectorConstruction* det)
: event_action(action), detector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  auto edep = step->GetTotalEnergyDeposit();
  if(edep <= 0.){
    return;
  }

  auto pre_step_point = step->GetPreStepPoint();
  if(!pre_step_point){
    return;
  }

  auto volume = pre_step_point->GetTouchableHandle()->GetVolume();
  if(!volume){
    return;
  }

  auto logical = volume->GetLogicalVolume();
  if(logical != detector->GetGAGGLogical()){
    return;
  }

  auto track = step->GetTrack();
  auto info = dynamic_cast<TrackInformation*>(track->GetUserInformation());
  if(!info){
    return;
  }

  if(event_action){
    event_action->AddEdep(info->GetPrimaryId(), edep);
  }
}
