#ifndef ActionInitialization_H
#define ActionInitialization_H 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"
#include "Rtypes.h"

class DetectorConstruction;
class RootIO;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(DetectorConstruction* detector, ULong64_t random_seed, const G4String& macro_file_name);
  ~ActionInitialization() override;

public:
  void Build() const override;

private:
  DetectorConstruction* detector = nullptr;
  ULong64_t random_seed = 0;
  G4String macro_file_name;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
