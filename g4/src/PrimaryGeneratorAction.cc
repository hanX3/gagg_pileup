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
#include "G4ios.hh"
#include "Randomize.hh"

#include "CLHEP/Random/RandPoisson.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct PrimarySample {
  G4double time = 0.0;
  G4double energy = 0.0;
  G4String particle_name = "gamma";
  G4String response_name = "other";
};

G4int SamplePoissonCount(G4double rate_hz, G4double time_window_s)
{
  const G4double mean = std::max(0.0, rate_hz)*time_window_s;
  return static_cast<G4int>(CLHEP::RandPoisson::shoot(mean));
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* action)
: event_action(action),
  use_disk_cone_source(DEFAULT_USE_DISK_CONE_SOURCE),
  aim_at_gagg(DEFAULT_AIM_AT_GAGG),
  source_disk_radius(SOURCE_DISK_RADIUS_MM*mm),
  source_cone_half_angle(SOURCE_EMISSION_CONE_HALF_ANGLE_DEG*deg),
  pileup_enabled(DEFAULT_PILEUP_ENABLED),
  brems_source_enabled(DEFAULT_BREMS_SOURCE_ENABLED),
  c12_capture_source_enabled(DEFAULT_C12_CAPTURE_SOURCE_ENABLED),
  alpha_p11b_source_enabled(DEFAULT_ALPHA_P11B_SOURCE_ENABLED),
  brems_source_rate_hz(DEFAULT_BREMS_SOURCE_RATE_HZ),
  c12_capture_source_rate_hz(DEFAULT_C12_CAPTURE_SOURCE_RATE_HZ),
  alpha_p11b_source_rate_hz(DEFAULT_ALPHA_P11B_SOURCE_RATE_HZ)
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

  pileup_messenger = new G4GenericMessenger(this, "/gagg/pileup/",
                                            "Poisson pileup source control.");

  auto& pileup_enable_cmd = pileup_messenger->DeclareProperty(
    "enable", pileup_enabled,
    "If true, one Geant4 event is one 1-ms waveform window with multiple primary photons."
  );
  pileup_enable_cmd.SetParameterName("enable", false);

  brems_source_messenger = new G4GenericMessenger(this, "/gagg/source/brems/",
                                                  "Thermal bremsstrahlung photon-source component.");

  auto& brems_enable_cmd = brems_source_messenger->DeclareProperty(
    "enable", brems_source_enabled,
    "Enable the Maxwellian thermal bremsstrahlung photon source. Fixed spectrum parameters are in Constants.hh."
  );
  brems_enable_cmd.SetParameterName("enable", false);

  auto& brems_rate_cmd = brems_source_messenger->DeclareProperty(
    "rateHz", brems_source_rate_hz,
    "Effective bremsstrahlung photon rate in Hz. Multiplicity is Poisson(rateHz * 1 ms)."
  );
  brems_rate_cmd.SetParameterName("rateHz", false);

  c12_capture_source_messenger = new G4GenericMessenger(this, "/gagg/source/c12Capture/",
                                                        "11B(p,gamma)12C capture/de-excitation photon-source component.");

  auto& c12_enable_cmd = c12_capture_source_messenger->DeclareProperty(
    "enable", c12_capture_source_enabled,
    "Enable the 12C capture/de-excitation photon source. The line-mixture model and line weights are fixed in Constants.hh."
  );
  c12_enable_cmd.SetParameterName("enable", false);

  auto& c12_rate_cmd = c12_capture_source_messenger->DeclareProperty(
    "rateHz", c12_capture_source_rate_hz,
    "Effective 12C capture/de-excitation photon rate in Hz. Multiplicity is Poisson(rateHz * 1 ms)."
  );
  c12_rate_cmd.SetParameterName("rateHz", false);


  alpha_p11b_source_messenger = new G4GenericMessenger(this, "/gagg/source/alphaP11B/",
                                                       "p-11B fusion alpha-particle source component.");

  auto& alpha_enable_cmd = alpha_p11b_source_messenger->DeclareProperty(
    "enable", alpha_p11b_source_enabled,
    "Enable the p-11B alpha source. The effective single-alpha energy spectrum is fixed in Constants.hh."
  );
  alpha_enable_cmd.SetParameterName("enable", false);

  auto& alpha_rate_cmd = alpha_p11b_source_messenger->DeclareProperty(
    "rateHz", alpha_p11b_source_rate_hz,
    "Effective p-11B alpha-particle rate in Hz. Multiplicity is Poisson(rateHz * 1 ms)."
  );
  alpha_rate_cmd.SetParameterName("rateHz", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete alpha_p11b_source_messenger;
  delete c12_capture_source_messenger;
  delete brems_source_messenger;
  delete pileup_messenger;
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
void PrimaryGeneratorAction::ConfigureSourcePositionAndDirection()
{
  if(use_disk_cone_source){
    const auto source_position = SampleDiskSourcePosition();
    particle_gun->SetParticlePosition(source_position);

    if(aim_at_gagg){
      particle_gun->SetParticleMomentumDirection(SampleAimAtGaggMomentumDirection(source_position));
    }else{
      particle_gun->SetParticleMomentumDirection(SampleConeMomentumDirection());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PrimaryGeneratorAction::SampleBremsstrahlungGammaEnergy() const
{
  // Parameterized p-11B plasma photon background:
  //   dN/dE proportional to E^{-1} exp(-E/kTe), Emin <= E <= Emax.
  // Fixed nominal parameters are in Constants.hh.
  //
  // Sampling method:
  //   Use log-uniform proposal q(E) proportional to 1/E, then accept with
  //   probability exp[-(E - Emin)/kTe].  This is much more efficient than
  //   uniform-energy rejection sampling for Emin << Emax.
  const G4double emin = std::max(1.0e-12*MeV, BREMS_EMIN);
  const G4double emax = std::max(emin, BREMS_EMAX);
  const G4double kte  = std::max(1.0e-12*MeV, BREMS_KTE);

  if(emax <= emin){
    return emin;
  }

  const G4double log_range = std::log(emax/emin);
  for(G4int iter = 0; iter < 10000; ++iter){
    const G4double e = emin*std::exp(log_range*G4UniformRand());
    const G4double accept = std::exp(-(e - emin)/kte);
    if(G4UniformRand() <= accept){
      return e;
    }
  }

  // Extremely unlikely fallback.  Return a valid energy inside the interval.
  return emin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PrimaryGeneratorAction::SampleC12CaptureGammaEnergy(G4String& response_name) const
{
  // Effective single-photon line-mixture approximation for 11B(p,gamma)12C.
  // This intentionally does not generate a true gamma1 + gamma_star cascade;
  // it samples one photon from the source component using photon-yield weights
  // derived from the gamma0/gamma1 branch ratio in Constants.hh.
  const G4double u = G4UniformRand();
  const G4double w0 = C12_CAPTURE_GAMMA0_PHOTON_WEIGHT;
  const G4double w1 = C12_CAPTURE_GAMMA1_PHOTON_WEIGHT;

  if(u < w0){
    response_name = "gamma_c12_capture_gamma0";
    return C12_CAPTURE_GAMMA0_ENERGY;
  }

  if(u < w0 + w1){
    response_name = "gamma_c12_capture_gamma1";
    return C12_CAPTURE_GAMMA1_ENERGY;
  }

  response_name = "gamma_c12_capture_gammaStar";
  return C12_CAPTURE_GAMMA_STAR_ENERGY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PrimaryGeneratorAction::SampleP11BAlphaEnergy(G4String& response_name) const
{
  // First-version effective single-alpha spectrum for
  //   p + 11B -> 3 alpha + 8.7 MeV.
  //
  // It is not a full three-body Dalitz generator.  The model treats the alpha
  // source as an effective particle-rate source:
  //   with probability 1/3: primary alpha, E = 3.76 MeV;
  //   with probability 2/3: secondary alpha from 8Be* decay, with
  //     E = Eboost + Estar + 2*sqrt(Eboost*Estar)*cos(theta),
  //     cos(theta) uniformly sampled in [-1, 1].
  if(G4UniformRand() < P11B_ALPHA_PRIMARY_FRACTION){
    response_name = "alpha_p11b_primary";
    return P11B_ALPHA_PRIMARY_ENERGY;
  }

  const G4double cos_theta = 2.0*G4UniformRand() - 1.0;
  const G4double energy = P11B_ALPHA_SECONDARY_EBOOST
                        + P11B_ALPHA_SECONDARY_ESTAR
                        + 2.0*std::sqrt(P11B_ALPHA_SECONDARY_EBOOST*P11B_ALPHA_SECONDARY_ESTAR)*cos_theta;

  response_name = "alpha_p11b_secondary";
  return std::max(0.0, energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::AttachPrimaryInformationToLastVertex(G4Event* event,
                                                                  G4int primary_id,
                                                                  const G4String& particle_name,
                                                                  const G4String& response_name) const
{
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if(!pileup_enabled){
    // Backward-compatible mode: one primary per Geant4 event.
    G4int primary_id = 0;

    auto particle_definition = particle_gun->GetParticleDefinition();
    G4String particle_name = "unknown";
    if(particle_definition){
      particle_name = particle_definition->GetParticleName();
    }
    G4String response_name = GetPrimaryResponseName(particle_name);

    ConfigureSourcePositionAndDirection();
    particle_gun->GeneratePrimaryVertex(event);
    AttachPrimaryInformationToLastVertex(event, primary_id, particle_name, response_name);
    return;
  }

  // Multi-component p-11B pileup mode.
  // Current components: bremsstrahlung photons, 12C capture/de-excitation
  // photons, and p-11B alpha particles.
  // One Geant4 event is one 1-ms waveform window.  Each source component has
  // its own effective particle/photon rate.  The generated primaries are sorted
  // by time, so primary_id is chronological across all source components.
  auto particle_table = G4ParticleTable::GetParticleTable();
  auto gamma = particle_table->FindParticle("gamma");
  auto alpha = particle_table->FindParticle("alpha");

  const G4double time_window_s = static_cast<G4double>(EVENT_TIME_LENGTH_PS)*1.0e-12;
  const G4double primary_time_window = static_cast<G4double>(EVENT_TIME_LENGTH_PS)*ps;

  std::vector<PrimarySample> primaries;

  if(brems_source_enabled){
    const G4int n_brems = SamplePoissonCount(brems_source_rate_hz, time_window_s);
    primaries.reserve(primaries.size() + n_brems);
    for(G4int i = 0; i < n_brems; ++i){
      PrimarySample sample;
      sample.time = primary_time_window*G4UniformRand();
      sample.energy = SampleBremsstrahlungGammaEnergy();
      sample.particle_name = "gamma";
      sample.response_name = "gamma_brems";
      primaries.push_back(sample);
    }
  }

  if(c12_capture_source_enabled){
    const G4int n_c12 = SamplePoissonCount(c12_capture_source_rate_hz, time_window_s);
    primaries.reserve(primaries.size() + n_c12);
    for(G4int i = 0; i < n_c12; ++i){
      PrimarySample sample;
      sample.time = primary_time_window*G4UniformRand();
      sample.energy = SampleC12CaptureGammaEnergy(sample.response_name);
      sample.particle_name = "gamma";
      primaries.push_back(sample);
    }
  }

  if(alpha_p11b_source_enabled){
    const G4int n_alpha = SamplePoissonCount(alpha_p11b_source_rate_hz, time_window_s);
    primaries.reserve(primaries.size() + n_alpha);
    for(G4int i = 0; i < n_alpha; ++i){
      PrimarySample sample;
      sample.time = primary_time_window*G4UniformRand();
      sample.energy = SampleP11BAlphaEnergy(sample.response_name);
      sample.particle_name = "alpha";
      primaries.push_back(sample);
    }
  }

  std::sort(primaries.begin(), primaries.end(),
            [](const PrimarySample& a, const PrimarySample& b){
              return a.time < b.time;
            });

  for(G4int primary_id = 0; primary_id < static_cast<G4int>(primaries.size()); ++primary_id){
    const auto& sample = primaries[primary_id];

    if(sample.particle_name == "alpha"){
      particle_gun->SetParticleDefinition(alpha);
    }else{
      particle_gun->SetParticleDefinition(gamma);
    }

    ConfigureSourcePositionAndDirection();
    particle_gun->SetParticleEnergy(sample.energy);
    particle_gun->SetParticleTime(sample.time);
    particle_gun->GeneratePrimaryVertex(event);
    AttachPrimaryInformationToLastVertex(event, primary_id, sample.particle_name, sample.response_name);
  }
}
