#ifndef PhysicsList_H
#define PhysicsList_H 1

#include "G4VModularPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList();
  ~PhysicsList() override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
