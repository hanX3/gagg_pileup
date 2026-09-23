#include "SiPMSD.hh"
#include "EventAction.hh"
#include "TrackInformation.hh"
#include "Constants.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiPMSD::SiPMSD(const G4String& name)
: G4VSensitiveDetector(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiPMSD::~SiPMSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto track = step->GetTrack();
  if(track->GetDefinition() != G4OpticalPhoton::Definition()){
    return false;
  }

  auto info = dynamic_cast<TrackInformation*>(track->GetUserInformation());
  if(!info){
    track->SetTrackStatus(fStopAndKill);
    return false;
  }

  auto pre_step_point = step->GetPreStepPoint();
  auto time = pre_step_point->GetGlobalTime();
  auto position = pre_step_point->GetPosition();

  auto event_action = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
  if(event_action){
    event_action->AddPhotonArrival(info->GetPrimaryId(), time, position);
  }

  // The first-version SiPM is a virtual absorbing sensitive layer.
  // Once an optical photon enters it, the arrival time is recorded and the photon is killed.
  track->SetTrackStatus(fStopAndKill);

  return true;
}
