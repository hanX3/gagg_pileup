#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
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

private:
  G4ParticleGun* particle_gun = nullptr;
  EventAction* event_action = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
