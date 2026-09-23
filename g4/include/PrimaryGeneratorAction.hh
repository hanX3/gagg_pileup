#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class EventAction;
struct RunInfoData;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(EventAction* event_action);
  ~PrimaryGeneratorAction() override;

public:
  void GeneratePrimaries(G4Event* event) override;
  void FillRunInfo(RunInfoData& data) const;

private:
  G4String GetPrimaryResponseName(const G4String& particle_name) const;
  G4ThreeVector SampleDiskSourcePosition() const;
  G4ThreeVector SampleConeMomentumDirection() const;
  G4ThreeVector SampleGaggTargetPosition() const;
  G4ThreeVector SampleAimAtGaggMomentumDirection(const G4ThreeVector& source_position) const;
  void ConfigureSourcePositionAndDirection();

  G4double SampleBremsstrahlungGammaEnergy() const;
  G4double SampleC12CaptureGammaEnergy(G4String& response_name) const;
  G4double SampleP11BAlphaEnergy(G4String& response_name) const;

  void AttachPrimaryInformationToLastVertex(G4Event* event,
                                            G4int primary_id,
                                            const G4String& particle_name,
                                            const G4String& response_name) const;

private:
  G4ParticleGun* particle_gun = nullptr;
  G4GenericMessenger* source_messenger = nullptr;
  G4GenericMessenger* pileup_messenger = nullptr;
  G4GenericMessenger* brems_source_messenger = nullptr;
  G4GenericMessenger* c12_capture_source_messenger = nullptr;
  G4GenericMessenger* alpha_p11b_source_messenger = nullptr;
  EventAction* event_action = nullptr;

  G4bool use_disk_cone_source = false;
  G4bool aim_at_gagg = false;
  G4double source_disk_radius = 0.0;
  G4double source_cone_half_angle = 0.0;

  // Pileup source:
  // one Geant4 event = one configurable waveform window of length T.
  // Each enabled p-11B source component is sampled independently
  // with Poisson(rate * T). All generated primaries are sorted by
  // emission time before assigning primary_id, so primary_id follows
  // chronological order inside the waveform window.
  G4bool pileup_enabled = false;

  G4bool brems_source_enabled = false;
  G4bool c12_capture_source_enabled = false;
  G4bool alpha_p11b_source_enabled = false;

  G4double brems_source_rate_hz = 0.0;
  G4double c12_capture_source_rate_hz = 0.0;
  G4double alpha_p11b_source_rate_hz = 0.0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
