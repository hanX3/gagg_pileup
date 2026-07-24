#ifndef PrimaryInformation_H
#define PrimaryInformation_H 1

#include "G4VUserPrimaryParticleInformation.hh"
#include "globals.hh"

#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PrimaryInformation : public G4VUserPrimaryParticleInformation
{
public:
  PrimaryInformation();
  PrimaryInformation(G4int primary_id,
                     const G4String& particle_name,
                     const G4String& response_name);
  ~PrimaryInformation() override;

public:
  void Print() const override;

public:
  void SetPrimaryId(G4int id) { primary_id = id; }
  void SetPrimaryParticleName(const G4String& name) { primary_particle_name = name; }
  void SetPrimaryResponseName(const G4String& name) { primary_response_name = name; }

  G4int GetPrimaryId() const { return primary_id; }
  const G4String& GetPrimaryParticleName() const { return primary_particle_name; }
  const G4String& GetPrimaryResponseName() const { return primary_response_name; }

private:
  G4int primary_id = -1;
  G4String primary_particle_name = "unknown";
  G4String primary_response_name = "other";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
