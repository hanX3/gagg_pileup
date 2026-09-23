#include "PhysicsList.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysicsList::PhysicsList()
{
  SetVerboseLevel(1);

  RegisterPhysics(new G4EmLivermorePhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4IonPhysics());
  RegisterPhysics(new G4OpticalPhysics());

  // Applies G4UserLimits/uStepMax to charged particles by default.
  // Do not call SetApplyToAll(true), so gamma rays and optical photons are not
  // forced into micron-scale zero-edep transport steps.
  RegisterPhysics(new G4StepLimiterPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysicsList::~PhysicsList()
{}
