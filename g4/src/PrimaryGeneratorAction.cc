#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "PrimaryInformation.hh"
#include "Constants.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* action)
: event_action(action)
{
  particle_gun = new G4ParticleGun(1);

  auto particle_table = G4ParticleTable::GetParticleTable();
  auto particle = particle_table->FindParticle("gamma");

  // Default source. These values can be overridden in the macro with /gun/*.
  particle_gun->SetParticleDefinition(particle);
  particle_gun->SetParticleEnergy(0.662*MeV);
  particle_gun->SetParticlePosition(G4ThreeVector(SOURCE_POSITION_X_MM*mm,
                                                  SOURCE_POSITION_Y_MM*mm,
                                                  SOURCE_POSITION_Z_MM*mm));
  particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particle_gun->SetParticleTime(0.*ps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particle_gun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String PrimaryGeneratorAction::GetPrimaryResponseName(const G4String& particle_name) const
{
  if(particle_name == "alpha"){
    return "alpha";
  }
  if(particle_name == "gamma" || particle_name == "e-" || particle_name == "e+"){
    return "gamma_electron";
  }

  return "other";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Current source model: one primary per Geant4 event.
  // Particle species, energy, position and direction are controlled by /gun/*
  // commands in the macro.
  G4int primary_id = 0;

  auto particle_definition = particle_gun->GetParticleDefinition();
  G4String particle_name = "unknown";
  if(particle_definition){
    particle_name = particle_definition->GetParticleName();
  }
  G4String response_name = GetPrimaryResponseName(particle_name);

  particle_gun->GeneratePrimaryVertex(event);

  auto vertex = event->GetPrimaryVertex(event->GetNumberOfPrimaryVertex()-1);
  if(vertex){
    auto primary_particle = vertex->GetPrimary();
    if(primary_particle){
      primary_particle->SetUserInformation(new PrimaryInformation(primary_id,
                                                                 particle_name,
                                                                 response_name));
    }
  }
}
