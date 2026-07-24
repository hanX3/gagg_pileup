#include "TrackInformation.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackInformation::TrackInformation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackInformation::TrackInformation(G4int id,
                                   const G4String& particle_name,
                                   const G4String& response_name)
: primary_id(id),
  primary_particle_name(particle_name),
  primary_response_name(response_name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackInformation::TrackInformation(const TrackInformation& right)
: G4VUserTrackInformation(right)
{
  primary_id = right.primary_id;
  primary_particle_name = right.primary_particle_name;
  primary_response_name = right.primary_response_name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackInformation::~TrackInformation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackInformation::Print() const
{
  G4cout << "TrackInformation: primary_id = " << primary_id
         << ", primary_particle_name = " << primary_particle_name
         << ", primary_response_name = " << primary_response_name << G4endl;
}
