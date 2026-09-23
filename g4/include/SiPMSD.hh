#ifndef SiPMSD_H
#define SiPMSD_H 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4TouchableHistory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SiPMSD : public G4VSensitiveDetector
{
public:
  SiPMSD(const G4String& name);
  ~SiPMSD() override;

public:
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
