#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;
class G4UserLimits;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

public:
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

public:
  G4LogicalVolume* GetGAGGLogical() const { return gagg_logical; }
  G4LogicalVolume* GetSiPMLogical() const { return sipm_logical; }

  void SetGaggMaxStep(G4double step);
  void SetMylarFrontMaxStep(G4double step);
  void SetMylarSideMaxStep(G4double step);

  G4double GetGaggMaxStep() const { return gagg_max_step; }
  G4double GetMylarFrontMaxStep() const { return mylar_front_max_step; }
  G4double GetMylarSideMaxStep() const { return mylar_side_max_step; }

private:
  void DefineMaterials();
  void ConstructOpticalProperties();
  void ConstructOpticalSurfaces();
  void PrintMaterialProperties();
  void ConstructStepLimitMessenger();
  void ApplyStepLimits();
  void PrintStepLimits();

private:
  G4Material* world_material = nullptr;
  G4Material* gagg_material = nullptr;
  G4Material* sipm_material = nullptr;
  G4Material* mylar_material = nullptr;

  G4GenericMessenger* step_messenger = nullptr;

  G4UserLimits* gagg_user_limits = nullptr;
  G4UserLimits* mylar_front_user_limits = nullptr;
  G4UserLimits* mylar_side_user_limits = nullptr;

  G4double gagg_max_step = 0.;
  G4double mylar_front_max_step = 0.;
  G4double mylar_side_max_step = 0.;

  G4LogicalVolume* world_logical = nullptr;
  G4LogicalVolume* gagg_logical = nullptr;
  G4LogicalVolume* sipm_logical = nullptr;
  G4LogicalVolume* mylar_front_logical = nullptr;
  G4LogicalVolume* mylar_xpos_logical = nullptr;
  G4LogicalVolume* mylar_xneg_logical = nullptr;
  G4LogicalVolume* mylar_ypos_logical = nullptr;
  G4LogicalVolume* mylar_yneg_logical = nullptr;

  G4VPhysicalVolume* gagg_physical = nullptr;
  G4VPhysicalVolume* sipm_physical = nullptr;
  G4VPhysicalVolume* mylar_front_physical = nullptr;
  G4VPhysicalVolume* mylar_xpos_physical = nullptr;
  G4VPhysicalVolume* mylar_xneg_physical = nullptr;
  G4VPhysicalVolume* mylar_ypos_physical = nullptr;
  G4VPhysicalVolume* mylar_yneg_physical = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
