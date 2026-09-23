#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class DetectorConstruction;
class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(EventAction* event_action, DetectorConstruction* detector);
  ~SteppingAction() override;

public:
  void UserSteppingAction(const G4Step* step) override;

private:
  EventAction* event_action = nullptr;
  DetectorConstruction* detector = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
