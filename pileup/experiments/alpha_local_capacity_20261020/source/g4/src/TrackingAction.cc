#include "TrackingAction.hh"
#include "TrackInformation.hh"
#include "PrimaryInformation.hh"
#include "EventAction.hh"

#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction(EventAction* action)
: event_action(action)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::~TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  if(track->GetUserInformation()){
    return;
  }

  if(track->GetParentID() == 0){
    G4int primary_id = track->GetTrackID() - 1;
    G4String particle_name = "unknown";
    G4String response_name = "other";

    auto primary_particle = track->GetDynamicParticle()->GetPrimaryParticle();
    if(primary_particle && primary_particle->GetUserInformation()){
      auto primary_info = dynamic_cast<PrimaryInformation*>(primary_particle->GetUserInformation());
      if(primary_info){
        primary_id = primary_info->GetPrimaryId();
        particle_name = primary_info->GetPrimaryParticleName();
        response_name = primary_info->GetPrimaryResponseName();
      }
    }

    const_cast<G4Track*>(track)->SetUserInformation(new TrackInformation(primary_id,
                                                                         particle_name,
                                                                         response_name));

    if(event_action){
      event_action->AddPrimary(primary_id,
                               particle_name,
                               response_name,
                               track->GetKineticEnergy(),
                               track->GetGlobalTime(),
                               track->GetPosition(),
                               track->GetMomentumDirection());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  auto info = dynamic_cast<TrackInformation*>(track->GetUserInformation());
  if(!info){
    return;
  }

  auto secondaries = fpTrackingManager->GimmeSecondaries();
  if(!secondaries){
    return;
  }

  for(auto secondary : *secondaries){
    if(!secondary){
      continue;
    }

    secondary->SetUserInformation(new TrackInformation(*info));

    if(secondary->GetDefinition() == G4OpticalPhoton::Definition()){
      auto creator = secondary->GetCreatorProcess();
      if(creator && creator->GetProcessName() == "Scintillation"){
        if(event_action){
          event_action->AddScintPhotonGenerated(info->GetPrimaryId());
        }
      }
    }
  }
}
