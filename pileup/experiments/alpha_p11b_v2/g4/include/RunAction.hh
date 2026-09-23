#ifndef RunAction_H
#define RunAction_H 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class RootIO;
class DetectorConstruction;
class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RunAction : public G4UserRunAction
{
public:
  RunAction(RootIO* root_io, DetectorConstruction* detector);
  ~RunAction() override;

public:
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

private:
  RootIO* root_io = nullptr;
  DetectorConstruction* detector = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
