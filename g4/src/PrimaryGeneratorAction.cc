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
#include "G4GenericMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* action)
: event_action(action),
  use_disk_cone_source(DEFAULT_USE_DISK_CONE_SOURCE),
  aim_at_gagg(DEFAULT_AIM_AT_GAGG),
  source_disk_radius(SOURCE_DISK_RADIUS_MM*mm),
  source_cone_half_angle(SOURCE_EMISSION_CONE_HALF_ANGLE_DEG*deg)
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

  source_messenger = new G4GenericMessenger(this, "/gagg/source/",
                                            "Primary source control for GAGG waveform simulations.");

  auto& use_cmd = source_messenger->DeclareProperty(
    "useDiskConeSource", use_disk_cone_source,
    "If true, override /gun/position and /gun/direction with a phi20-mm disk source and a limited emission cone around +z."
  );
  use_cmd.SetParameterName("useDiskConeSource", false);

  auto& aim_cmd = source_messenger->DeclareProperty(
    "aimAtGagg", aim_at_gagg,
    "If true and useDiskConeSource is true, sample a random target point on the GAGG front face and aim the primary at it."
  );
  aim_cmd.SetParameterName("aimAtGagg", false);

  auto& radius_cmd = source_messenger->DeclarePropertyWithUnit(
    "diskRadius", "mm", source_disk_radius,
    "Disk-source radius. Default is SOURCE_DISK_RADIUS_MM in Constants.hh."
  );
  radius_cmd.SetParameterName("diskRadius", false);

  auto& angle_cmd = source_messenger->DeclarePropertyWithUnit(
    "coneHalfAngle", "deg", source_cone_half_angle,
    "Emission-cone half angle around +z. Default is SOURCE_EMISSION_CONE_HALF_ANGLE_DEG in Constants.hh."
  );
  angle_cmd.SetParameterName("coneHalfAngle", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete source_messenger;
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
G4ThreeVector PrimaryGeneratorAction::SampleDiskSourcePosition() const
{
  const G4double radius = std::max(0.0, source_disk_radius);
  const G4double r = radius*std::sqrt(G4UniformRand());
  const G4double phi = twopi*G4UniformRand();

  return G4ThreeVector(SOURCE_POSITION_X_MM*mm + r*std::cos(phi),
                       SOURCE_POSITION_Y_MM*mm + r*std::sin(phi),
                       SOURCE_POSITION_Z_MM*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector PrimaryGeneratorAction::SampleConeMomentumDirection() const
{
  const G4double theta_max = std::max(0.0, std::min(pi, source_cone_half_angle));
  if(theta_max <= 0.0){
    return G4ThreeVector(0., 0., 1.);
  }

  // Uniform in solid angle inside the cone.
  const G4double cos_theta_max = std::cos(theta_max);
  const G4double cos_theta = cos_theta_max + (1.0 - cos_theta_max)*G4UniformRand();
  const G4double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta*cos_theta));
  const G4double phi = twopi*G4UniformRand();

  return G4ThreeVector(sin_theta*std::cos(phi),
                       sin_theta*std::sin(phi),
                       cos_theta).unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector PrimaryGeneratorAction::SampleGaggTargetPosition() const
{
  const G4double x = SOURCE_TARGET_CENTER_X_MM*mm
                   + (G4UniformRand() - 0.5)*SOURCE_TARGET_SIZE_X_MM*mm;
  const G4double y = SOURCE_TARGET_CENTER_Y_MM*mm
                   + (G4UniformRand() - 0.5)*SOURCE_TARGET_SIZE_Y_MM*mm;
  const G4double z = SOURCE_TARGET_Z_MM*mm;

  return G4ThreeVector(x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector PrimaryGeneratorAction::SampleAimAtGaggMomentumDirection(const G4ThreeVector& source_position) const
{
  const auto target_position = SampleGaggTargetPosition();
  auto direction = target_position - source_position;

  if(direction.mag2() <= 0.0){
    return G4ThreeVector(0., 0., 1.);
  }

  return direction.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Current source model: one primary per Geant4 event.
  // Particle species and energy are controlled by /gun/* commands in the macro.
  // If /gagg/source/useDiskConeSource is true, the source position and momentum
  // direction are sampled here event by event.
  G4int primary_id = 0;

  auto particle_definition = particle_gun->GetParticleDefinition();
  G4String particle_name = "unknown";
  if(particle_definition){
    particle_name = particle_definition->GetParticleName();
  }
  G4String response_name = GetPrimaryResponseName(particle_name);

  if(use_disk_cone_source){
    const auto source_position = SampleDiskSourcePosition();
    particle_gun->SetParticlePosition(source_position);

    if(aim_at_gagg){
      particle_gun->SetParticleMomentumDirection(SampleAimAtGaggMomentumDirection(source_position));
    }else{
      particle_gun->SetParticleMomentumDirection(SampleConeMomentumDirection());
    }
  }

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
