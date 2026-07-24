#ifndef EventAction_H
#define EventAction_H 1

#include "DataStructure.hh"
#include "Constants.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class RootIO;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class EventAction : public G4UserEventAction
{
public:
  EventAction(RootIO* root_io);
  ~EventAction() override;

public:
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

public:
  void AddPrimary(G4int primary_id,
                  const G4String& particle_name,
                  const G4String& response_name,
                  G4double energy,
                  G4double time,
                  const G4ThreeVector& position,
                  const G4ThreeVector& direction);

  void AddEdep(G4int primary_id, G4double edep);
  void AddScintPhotonGenerated(G4int primary_id);
  void AddPhotonArrival(G4int primary_id,
                        G4double time,
                        const G4ThreeVector& position);

  G4int GetNPrimary() const { return primary_photon_data_vec.size(); }

private:
  G4bool IsValidPrimaryId(G4int primary_id) const;

private:
  RootIO* root_io = nullptr;

  WaveformEventData waveform_event_data;
  std::vector<PrimaryPhotonData> primary_photon_data_vec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
