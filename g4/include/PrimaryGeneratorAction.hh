#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(EventAction* event_action);
  ~PrimaryGeneratorAction() override;

public:
  void GeneratePrimaries(G4Event* event) override;

private:
  G4String GetPrimaryResponseName(const G4String& particle_name) const;
  G4ThreeVector SampleDiskSourcePosition() const;
  G4ThreeVector SampleConeMomentumDirection() const;
  G4ThreeVector SampleGaggTargetPosition() const;
  G4ThreeVector SampleAimAtGaggMomentumDirection(const G4ThreeVector& source_position) const;

private:
  G4ParticleGun* particle_gun = nullptr;
  G4GenericMessenger* source_messenger = nullptr;
  EventAction* event_action = nullptr;

  G4bool use_disk_cone_source = false;
  G4bool aim_at_gagg = false;
  G4double source_disk_radius = 0.0;
  G4double source_cone_half_angle = 0.0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
