#include "PrimaryInformation.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryInformation::PrimaryInformation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryInformation::PrimaryInformation(G4int id,
                                       const G4String& particle_name,
                                       const G4String& response_name)
: primary_id(id),
  primary_particle_name(particle_name),
  primary_response_name(response_name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryInformation::~PrimaryInformation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryInformation::Print() const
{
  G4cout << "PrimaryInformation: primary_id = " << primary_id
         << ", primary_particle_name = " << primary_particle_name
         << ", primary_response_name = " << primary_response_name << G4endl;
}
