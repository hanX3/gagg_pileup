#include "DetectorConstruction.hh"
#include "Constants.hh"
#include "SiPMSD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4IonisParamMat.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
  gagg_max_step = DEFAULT_GAGG_MAX_STEP_UM*um;
  mylar_front_max_step = DEFAULT_MYLAR_FRONT_MAX_STEP_UM*um;
  mylar_side_max_step = DEFAULT_MYLAR_SIDE_MAX_STEP_UM*um;

  ConstructStepLimitMessenger();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  delete step_messenger;
  delete gagg_user_limits;
  delete mylar_front_user_limits;
  delete mylar_side_user_limits;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructStepLimitMessenger()
{
  step_messenger = new G4GenericMessenger(this, "/gagg/step/", "GAGG step-limit control");

  auto& gagg_cmd = step_messenger->DeclareMethodWithUnit(
    "setGaggMaxStep", "um",
    &DetectorConstruction::SetGaggMaxStep,
    "Set maximum charged-particle step length in the GAGG crystal."
  );
  gagg_cmd.SetParameterName("gaggMaxStep", false);
  gagg_cmd.SetRange("gaggMaxStep>0.");

  auto& front_cmd = step_messenger->DeclareMethodWithUnit(
    "setMylarFrontMaxStep", "um",
    &DetectorConstruction::SetMylarFrontMaxStep,
    "Set maximum charged-particle step length in the front Mylar window."
  );
  front_cmd.SetParameterName("mylarFrontMaxStep", false);
  front_cmd.SetRange("mylarFrontMaxStep>0.");

  auto& side_cmd = step_messenger->DeclareMethodWithUnit(
    "setMylarSideMaxStep", "um",
    &DetectorConstruction::SetMylarSideMaxStep,
    "Set maximum charged-particle step length in the side Mylar films."
  );
  side_cmd.SetParameterName("mylarSideMaxStep", false);
  side_cmd.SetRange("mylarSideMaxStep>0.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetGaggMaxStep(G4double step)
{
  gagg_max_step = step;
  ApplyStepLimits();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMylarFrontMaxStep(G4double step)
{
  mylar_front_max_step = step;
  ApplyStepLimits();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMylarSideMaxStep(G4double step)
{
  mylar_side_max_step = step;
  ApplyStepLimits();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ApplyStepLimits()
{
  if(gagg_logical){
    if(!gagg_user_limits){
      gagg_user_limits = new G4UserLimits(gagg_max_step);
    }else{
      gagg_user_limits->SetMaxAllowedStep(gagg_max_step);
    }
    gagg_logical->SetUserLimits(gagg_user_limits);
  }

  if(mylar_front_logical){
    if(!mylar_front_user_limits){
      mylar_front_user_limits = new G4UserLimits(mylar_front_max_step);
    }else{
      mylar_front_user_limits->SetMaxAllowedStep(mylar_front_max_step);
    }
    mylar_front_logical->SetUserLimits(mylar_front_user_limits);
  }

  if(mylar_xpos_logical || mylar_xneg_logical || mylar_ypos_logical || mylar_yneg_logical){
    if(!mylar_side_user_limits){
      mylar_side_user_limits = new G4UserLimits(mylar_side_max_step);
    }else{
      mylar_side_user_limits->SetMaxAllowedStep(mylar_side_max_step);
    }

    if(mylar_xpos_logical){ mylar_xpos_logical->SetUserLimits(mylar_side_user_limits); }
    if(mylar_xneg_logical){ mylar_xneg_logical->SetUserLimits(mylar_side_user_limits); }
    if(mylar_ypos_logical){ mylar_ypos_logical->SetUserLimits(mylar_side_user_limits); }
    if(mylar_yneg_logical){ mylar_yneg_logical->SetUserLimits(mylar_side_user_limits); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintStepLimits()
{
  G4cout << G4endl;
  G4cout << "========== Charged-Particle Max Step Check ==========" << G4endl;
  G4cout << "GAGG max step        : " << gagg_max_step/um << " um" << G4endl;
  G4cout << "Front Mylar max step : " << mylar_front_max_step/um << " um" << G4endl;
  G4cout << "Side Mylar max step  : " << mylar_side_max_step/um << " um" << G4endl;
  G4cout << "Note: with the default G4StepLimiterPhysics configuration, " << G4endl;
  G4cout << "      these limits apply to charged particles, not to gamma rays or optical photons." << G4endl;
  G4cout << "=====================================================" << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials()
{
  auto nist = G4NistManager::Instance();

  world_material = nist->FindOrBuildMaterial("G4_Galactic");
  mylar_material = nist->FindOrBuildMaterial("G4_MYLAR");

  auto elGd = nist->FindOrBuildElement("Gd");
  auto elAl = nist->FindOrBuildElement("Al");
  auto elGa = nist->FindOrBuildElement("Ga");
  auto elO  = nist->FindOrBuildElement("O");

  gagg_material = new G4Material("GAGG", GAGG_DENSITY_G_CM3*g/cm3, 4);
  gagg_material->AddElement(elGd, 3);
  gagg_material->AddElement(elAl, 2);
  gagg_material->AddElement(elGa, 3);
  gagg_material->AddElement(elO,  12);

  auto elSi = nist->FindOrBuildElement("Si");
  sipm_material = new G4Material("VirtualSiPM", 2.33*g/cm3, 1);
  sipm_material->AddElement(elSi, 1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructOpticalProperties()
{
  // First-version simplified optical properties.
  // The GAGG emission component is centered roughly around 530 nm.
  const G4int n_entries = 8;
  G4double photon_energy[n_entries] = {
    2.00*eV, 2.15*eV, 2.25*eV, 2.34*eV,
    2.43*eV, 2.55*eV, 2.75*eV, 3.00*eV
  };

  G4double rindex_gagg[n_entries] = {
    GAGG_RINDEX, GAGG_RINDEX, GAGG_RINDEX, GAGG_RINDEX,
    GAGG_RINDEX, GAGG_RINDEX, GAGG_RINDEX, GAGG_RINDEX
  };

  G4double abs_length_gagg[n_entries] = {
    GAGG_ABS_LENGTH_CM*cm, GAGG_ABS_LENGTH_CM*cm,
    GAGG_ABS_LENGTH_CM*cm, GAGG_ABS_LENGTH_CM*cm,
    GAGG_ABS_LENGTH_CM*cm, GAGG_ABS_LENGTH_CM*cm,
    GAGG_ABS_LENGTH_CM*cm, GAGG_ABS_LENGTH_CM*cm
  };

  G4double scint_fast[n_entries] = {
    0.05, 0.20, 0.55, 1.00,
    0.55, 0.20, 0.05, 0.01
  };

  G4double scint_slow[n_entries] = {
    0.05, 0.20, 0.55, 1.00,
    0.55, 0.20, 0.05, 0.01
  };

  auto gagg_mpt = new G4MaterialPropertiesTable();
  gagg_mpt->AddProperty("RINDEX", photon_energy, rindex_gagg, n_entries);
  gagg_mpt->AddProperty("ABSLENGTH", photon_energy, abs_length_gagg, n_entries);
  gagg_mpt->AddProperty("SCINTILLATIONCOMPONENT1", photon_energy, scint_fast, n_entries);
  gagg_mpt->AddProperty("SCINTILLATIONCOMPONENT2", photon_energy, scint_slow, n_entries);

  gagg_mpt->AddConstProperty("SCINTILLATIONYIELD", GAGG_SCINT_YIELD_PER_MEV/MeV);
  gagg_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  gagg_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", GAGG_FAST_TIME_NS*ns);
  gagg_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", GAGG_SLOW_TIME_NS*ns);
  gagg_mpt->AddConstProperty("SCINTILLATIONYIELD1", GAGG_FAST_FRACTION);
  gagg_mpt->AddConstProperty("SCINTILLATIONYIELD2", GAGG_SLOW_FRACTION);

  gagg_material->SetMaterialPropertiesTable(gagg_mpt);

  // Enable Birks-type scintillation quenching in Geant4 optical photon generation.
  // In the current ROOT scheme, edep_total_MeV is the truth energy and the
  // optical-photon branches carry the detector optical response.
  if(gagg_material->GetIonisation()){
    gagg_material->GetIonisation()->SetBirksConstant(GAGG_BIRKS_CONSTANT_MM_PER_MEV*mm/MeV);
  }

  // Virtual SiPM uses the same refractive index as GAGG to suppress Fresnel reflection
  // at the GAGG/SiPM interface in the first-version model.
  G4double rindex_sipm[n_entries] = {
    VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX,
    VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX, VIRTUAL_SIPM_RINDEX
  };
  auto sipm_mpt = new G4MaterialPropertiesTable();
  sipm_mpt->AddProperty("RINDEX", photon_energy, rindex_sipm, n_entries);
  sipm_material->SetMaterialPropertiesTable(sipm_mpt);

  G4double rindex_vacuum[n_entries] = {
    VACUUM_RINDEX, VACUUM_RINDEX, VACUUM_RINDEX, VACUUM_RINDEX,
    VACUUM_RINDEX, VACUUM_RINDEX, VACUUM_RINDEX, VACUUM_RINDEX
  };
  auto vacuum_mpt = new G4MaterialPropertiesTable();
  vacuum_mpt->AddProperty("RINDEX", photon_energy, rindex_vacuum, n_entries);
  world_material->SetMaterialPropertiesTable(vacuum_mpt);

  G4double rindex_mylar[n_entries] = {
    MYLAR_RINDEX, MYLAR_RINDEX, MYLAR_RINDEX, MYLAR_RINDEX,
    MYLAR_RINDEX, MYLAR_RINDEX, MYLAR_RINDEX, MYLAR_RINDEX
  };

  G4double abs_length_mylar[n_entries] = {
    MYLAR_ABS_LENGTH_CM*cm, MYLAR_ABS_LENGTH_CM*cm,
    MYLAR_ABS_LENGTH_CM*cm, MYLAR_ABS_LENGTH_CM*cm,
    MYLAR_ABS_LENGTH_CM*cm, MYLAR_ABS_LENGTH_CM*cm,
    MYLAR_ABS_LENGTH_CM*cm, MYLAR_ABS_LENGTH_CM*cm
  };

  auto mylar_mpt = new G4MaterialPropertiesTable();
  mylar_mpt->AddProperty("RINDEX", photon_energy, rindex_mylar, n_entries);
  mylar_mpt->AddProperty("ABSLENGTH", photon_energy, abs_length_mylar, n_entries);
  mylar_material->SetMaterialPropertiesTable(mylar_mpt);

  PrintMaterialProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintMaterialProperties()
{
  G4cout << G4endl;
  G4cout << "========== Material and Optical Properties Check ==========" << G4endl;

  G4Material* materials[4] = {
    world_material,
    gagg_material,
    mylar_material,
    sipm_material
  };

  for(G4int i=0;i<4;i++){
    auto material = materials[i];

    if(!material){
      G4cout << G4endl;
      G4cout << "Material pointer is null." << G4endl;
      continue;
    }

    G4cout << G4endl;
    G4cout << "Material: " << material->GetName() << G4endl;
    G4cout << "Density : " << material->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
    if(material->GetIonisation()){
      G4cout << "Birks constant : "
             << material->GetIonisation()->GetBirksConstant()/(mm/MeV)
             << " mm/MeV" << G4endl;
    }

    auto mpt = material->GetMaterialPropertiesTable();
    if(!mpt){
      G4cout << "Optical MPT: none" << G4endl;
      continue;
    }

    G4cout << "Optical MPT exists. DumpTable():" << G4endl;
    mpt->DumpTable();
  }

  G4cout << "===========================================================" << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructOpticalSurfaces()
{
  // Side Mylar films and the front Mylar window are real dielectric layers.
  // These optical surfaces are not physical objects by themselves; they only
  // define the boundary treatment at selected GAGG/Mylar interfaces.
  // No artificial high reflectivity is imposed for pure Mylar.
  auto side_mylar_surface = new G4OpticalSurface("GAGGSideMylarSurface");
  side_mylar_surface->SetType(dielectric_dielectric);
  side_mylar_surface->SetModel(unified);
  side_mylar_surface->SetFinish(polished);

  if(gagg_physical && mylar_xpos_physical){
    new G4LogicalBorderSurface("GAGG_Mylar_XPos_Surface", gagg_physical, mylar_xpos_physical, side_mylar_surface);
  }
  if(gagg_physical && mylar_xneg_physical){
    new G4LogicalBorderSurface("GAGG_Mylar_XNeg_Surface", gagg_physical, mylar_xneg_physical, side_mylar_surface);
  }
  if(gagg_physical && mylar_ypos_physical){
    new G4LogicalBorderSurface("GAGG_Mylar_YPos_Surface", gagg_physical, mylar_ypos_physical, side_mylar_surface);
  }
  if(gagg_physical && mylar_yneg_physical){
    new G4LogicalBorderSurface("GAGG_Mylar_YNeg_Surface", gagg_physical, mylar_yneg_physical, side_mylar_surface);
  }

  auto front_surface = new G4OpticalSurface("GAGGMylarFrontSurface");
  front_surface->SetType(dielectric_dielectric);
  front_surface->SetModel(unified);
  front_surface->SetFinish(polished);

  if(gagg_physical && mylar_front_physical){
    new G4LogicalBorderSurface("GAGG_MylarFront_Surface", gagg_physical, mylar_front_physical, front_surface);
  }

  G4cout << G4endl;
  G4cout << "========== Optical Surface Check ==========" << G4endl;
  G4cout << "GAGG -> side Mylar surface:" << G4endl;
  side_mylar_surface->DumpInfo();
  G4cout << "GAGG -> front Mylar surface:" << G4endl;
  front_surface->DumpInfo();
  G4cout << "===========================================" << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  ConstructOpticalProperties();

  auto world_solid = new G4Box("World",
                               0.5*WORLD_SIZE_MM*mm,
                               0.5*WORLD_SIZE_MM*mm,
                               0.5*WORLD_SIZE_MM*mm);
  world_logical = new G4LogicalVolume(world_solid, world_material, "World_logical");
  auto world_physical = new G4PVPlacement(nullptr, G4ThreeVector(), world_logical,
                                          "World", nullptr, false, 0, true);

  auto gagg_solid = new G4Box("GAGG",
                              0.5*GAGG_SIZE_X_MM*mm,
                              0.5*GAGG_SIZE_Y_MM*mm,
                              0.5*GAGG_SIZE_Z_MM*mm);
  gagg_logical = new G4LogicalVolume(gagg_solid, gagg_material, "GAGG_logical");
  gagg_physical = new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,GAGG_CENTER_Z_MM*mm), gagg_logical,
                                    "GAGG", world_logical, false, 0, true);

  auto sipm_solid = new G4Box("SiPM",
                              0.5*SIPM_SIZE_X_MM*mm,
                              0.5*SIPM_SIZE_Y_MM*mm,
                              0.5*SIPM_SIZE_Z_MM*mm);
  sipm_logical = new G4LogicalVolume(sipm_solid, sipm_material, "SiPM_logical");
  sipm_physical = new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,SIPM_CENTER_Z_MM*mm), sipm_logical,
                                    "SiPM", world_logical, false, 0, true);

  auto mylar_front_solid = new G4Box("MylarFront",
                                     0.5*MYLAR_FRONT_SIZE_X_MM*mm,
                                     0.5*MYLAR_FRONT_SIZE_Y_MM*mm,
                                     0.5*MYLAR_FRONT_SIZE_Z_MM*mm);
  mylar_front_logical = new G4LogicalVolume(mylar_front_solid, mylar_material, "MylarFront_logical");
  mylar_front_physical = new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,MYLAR_FRONT_CENTER_Z_MM*mm), mylar_front_logical,
                                           "MylarFront", world_logical, false, 0, true);

  auto mylar_x_solid = new G4Box("Mylar_XSide",
                                 0.5*MYLAR_SIDE_THICKNESS_MM*mm,
                                 0.5*GAGG_SIZE_Y_MM*mm,
                                 0.5*GAGG_SIZE_Z_MM*mm);
  mylar_xpos_logical = new G4LogicalVolume(mylar_x_solid, mylar_material, "Mylar_XPos_logical");
  mylar_xneg_logical = new G4LogicalVolume(mylar_x_solid, mylar_material, "Mylar_XNeg_logical");

  mylar_xpos_physical = new G4PVPlacement(nullptr,
                                          G4ThreeVector((0.5*GAGG_SIZE_X_MM + 0.5*MYLAR_SIDE_THICKNESS_MM)*mm,
                                                        0.,
                                                        GAGG_CENTER_Z_MM*mm),
                                          mylar_xpos_logical,
                                          "Mylar_XPos", world_logical, false, 0, true);

  mylar_xneg_physical = new G4PVPlacement(nullptr,
                                          G4ThreeVector((-0.5*GAGG_SIZE_X_MM - 0.5*MYLAR_SIDE_THICKNESS_MM)*mm,
                                                        0.,
                                                        GAGG_CENTER_Z_MM*mm),
                                          mylar_xneg_logical,
                                          "Mylar_XNeg", world_logical, false, 0, true);

  auto mylar_y_solid = new G4Box("Mylar_YSide",
                                 0.5*(GAGG_SIZE_X_MM + 2.0*MYLAR_SIDE_THICKNESS_MM)*mm,
                                 0.5*MYLAR_SIDE_THICKNESS_MM*mm,
                                 0.5*GAGG_SIZE_Z_MM*mm);
  mylar_ypos_logical = new G4LogicalVolume(mylar_y_solid, mylar_material, "Mylar_YPos_logical");
  mylar_yneg_logical = new G4LogicalVolume(mylar_y_solid, mylar_material, "Mylar_YNeg_logical");

  mylar_ypos_physical = new G4PVPlacement(nullptr,
                                          G4ThreeVector(0.,
                                                        (0.5*GAGG_SIZE_Y_MM + 0.5*MYLAR_SIDE_THICKNESS_MM)*mm,
                                                        GAGG_CENTER_Z_MM*mm),
                                          mylar_ypos_logical,
                                          "Mylar_YPos", world_logical, false, 0, true);

  mylar_yneg_physical = new G4PVPlacement(nullptr,
                                          G4ThreeVector(0.,
                                                        (-0.5*GAGG_SIZE_Y_MM - 0.5*MYLAR_SIDE_THICKNESS_MM)*mm,
                                                        GAGG_CENTER_Z_MM*mm),
                                          mylar_yneg_logical,
                                          "Mylar_YNeg", world_logical, false, 0, true);

  ConstructOpticalSurfaces();
  ApplyStepLimits();
  PrintStepLimits();

  auto world_vis = new G4VisAttributes();
  world_vis->SetVisibility(false);
  world_logical->SetVisAttributes(world_vis);

  auto gagg_vis = new G4VisAttributes(G4Colour(0.1,0.8,0.1,0.35));
  gagg_vis->SetForceSolid(true);
  gagg_logical->SetVisAttributes(gagg_vis);

  auto sipm_vis = new G4VisAttributes(G4Colour(0.1,0.1,0.9,0.50));
  sipm_vis->SetForceSolid(true);
  sipm_logical->SetVisAttributes(sipm_vis);

  auto mylar_vis = new G4VisAttributes(G4Colour(0.9,0.9,0.1,0.35));
  mylar_vis->SetForceSolid(true);
  mylar_front_logical->SetVisAttributes(mylar_vis);

  auto mylar_side_vis = new G4VisAttributes(G4Colour(0.9,0.9,0.1,0.25));
  mylar_side_vis->SetForceSolid(true);
  mylar_xpos_logical->SetVisAttributes(mylar_side_vis);
  mylar_xneg_logical->SetVisAttributes(mylar_side_vis);
  mylar_ypos_logical->SetVisAttributes(mylar_side_vis);
  mylar_yneg_logical->SetVisAttributes(mylar_side_vis);

  return world_physical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  auto sd_manager = G4SDManager::GetSDMpointer();

  auto sipm_sd = new SiPMSD("SiPMSD");
  sd_manager->AddNewDetector(sipm_sd);

  if(sipm_logical){
    sipm_logical->SetSensitiveDetector(sipm_sd);
  }
}
