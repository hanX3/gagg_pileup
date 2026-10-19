#ifndef TrackingAction_H
#define TrackingAction_H 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class EventAction;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class TrackingAction : public G4UserTrackingAction
{
public:
  TrackingAction(EventAction* event_action);
  ~TrackingAction() override;

public:
  void PreUserTrackingAction(const G4Track* track) override;
  void PostUserTrackingAction(const G4Track* track) override;

private:
  EventAction* event_action = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
