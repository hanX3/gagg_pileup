#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"

#include "G4EmCalculator.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
// Settings
// -----------------------------------------------------------------------------
namespace Settings {
    constexpr double kDensityGcm3 = 6.63;
    constexpr double kEminMeV = 0.01;
    constexpr double kEmaxMeV = 10.0;
    constexpr int    kNPoints = 500;

    const std::string kElectronOut = "../electron_gagg_stopping.csv";
    const std::string kAlphaOut    = "../alpha_gagg_stopping.csv";
    const std::string kO16Out      = "../o16_gagg_stopping.csv";
}


G4Material* BuildGAGG()
{
    auto nist = G4NistManager::Instance();

    auto Gd = nist->FindOrBuildElement("Gd");
    auto Al = nist->FindOrBuildElement("Al");
    auto Ga = nist->FindOrBuildElement("Ga");
    auto O  = nist->FindOrBuildElement("O");

    auto gagg = new G4Material("GAGG", Settings::kDensityGcm3 * g/cm3, 4);
    gagg->AddElement(Gd, 3);
    gagg->AddElement(Al, 2);
    gagg->AddElement(Ga, 3);
    gagg->AddElement(O,  12);

    return gagg;
}


class DetectorConstruction final : public G4VUserDetectorConstruction
{
public:
    G4VPhysicalVolume* Construct() override
    {
        auto gagg = BuildGAGG();

        auto solidWorld = new G4Box("World", 10.0*cm, 10.0*cm, 10.0*cm);
        auto logicWorld = new G4LogicalVolume(solidWorld, gagg, "WorldLV");

        return new G4PVPlacement(
            nullptr,
            G4ThreeVector(),
            logicWorld,
            "WorldPV",
            nullptr,
            false,
            0,
            true
        );
    }
};


class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction()
    {
        fGun = new G4ParticleGun(1);

        auto particle = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
        if (!particle) {
            throw std::runtime_error("Cannot find geantino.");
        }

        fGun->SetParticleDefinition(particle);
        fGun->SetParticleEnergy(1.0*MeV);
        fGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
        fGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
    }

    ~PrimaryGeneratorAction() override
    {
        delete fGun;
    }

    void GeneratePrimaries(G4Event* event) override
    {
        fGun->GeneratePrimaryVertex(event);
    }

private:
    G4ParticleGun* fGun = nullptr;
};


std::vector<G4double> MakeLogEnergyGrid(double eMinMeV,
                                        double eMaxMeV,
                                        int nPoints)
{
    std::vector<G4double> energies;
    energies.reserve(nPoints);

    const double logMin = std::log10(eMinMeV);
    const double logMax = std::log10(eMaxMeV);

    for (int i = 0; i < nPoints; ++i) {
        const double u = static_cast<double>(i) / static_cast<double>(nPoints - 1);
        const double eMeV = std::pow(10.0, logMin + u * (logMax - logMin));
        energies.push_back(eMeV * MeV);
    }

    return energies;
}


void PrintMaterialCutsCouples()
{
    auto table = G4ProductionCutsTable::GetProductionCutsTable();

    std::cout << "Material-cuts couples:" << std::endl;

    const std::size_t n = table->GetTableSize();
    for (std::size_t i = 0; i < n; ++i) {
        auto couple = table->GetMaterialCutsCouple(i);
        if (!couple) continue;

        auto mat = couple->GetMaterial();
        if (!mat) continue;

        std::cout << "  " << i
                  << " : " << mat->GetName()
                  << "  density = " << mat->GetDensity()/(g/cm3)
                  << " g/cm3"
                  << std::endl;
    }
}


void WriteStoppingTableForParticle(const std::string& filename,
                                   G4ParticleDefinition* particle,
                                   const std::string& particleLabel,
                                   const std::string& processName,
                                   G4Material* material,
                                   const std::vector<G4double>& energies)
{
    if (!particle) {
        throw std::runtime_error("Null particle definition for: " + particleLabel);
    }

    G4EmCalculator emCal;

    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    fout << "E_MeV,S_mass_MeV_cm2_g\n";

    std::cout << "Writing " << filename
              << " using TestEm0-style G4EmCalculator::ComputeDEDX"
              << "  particle=" << particleLabel
              << "  process=" << processName
              << std::endl;

    for (const auto E : energies) {
        // TestEm0-style process-specific dE/dx query.
        // For electron: process = eIoni
        // For alpha/ions: process = ionIoni
        // E is the total kinetic energy of the particle, not MeV/u.
        const G4double cut = DBL_MAX;
        const G4double dedx = emCal.ComputeDEDX(
            E,
            particle,
            processName,
            material,
            cut
        );

        const G4double sMass = dedx / material->GetDensity();
        const G4double sMassOut = sMass / (MeV * cm2 / g);

        fout << std::setprecision(12)
             << E / MeV << ","
             << sMassOut << "\n";
    }

    fout.close();
    std::cout << "Saved: " << filename << std::endl;
}


void WriteStoppingTable(const std::string& filename,
                        const std::string& particleName,
                        const std::string& processName,
                        G4Material* material,
                        const std::vector<G4double>& energies)
{
    auto particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    if (!particle) {
        throw std::runtime_error("Cannot find particle: " + particleName);
    }

    WriteStoppingTableForParticle(
        filename,
        particle,
        particleName,
        processName,
        material,
        energies
    );
}


G4ParticleDefinition* GetIonDefinition(G4int Z, G4int A, G4double excitationEnergy = 0.0)
{
    auto particleTable = G4ParticleTable::GetParticleTable();
    if (!particleTable) {
        throw std::runtime_error("Cannot access G4ParticleTable.");
    }

    auto ionTable = particleTable->GetIonTable();
    if (!ionTable) {
        throw std::runtime_error("Cannot access G4IonTable.");
    }

    auto ion = ionTable->GetIon(Z, A, excitationEnergy);
    if (!ion) {
        throw std::runtime_error(
            "Cannot create ion definition for Z=" + std::to_string(Z) +
            ", A=" + std::to_string(A)
        );
    }

    return ion;
}


class RunAction final : public G4UserRunAction
{
public:
    void BeginOfRunAction(const G4Run*) override
    {
        PrintMaterialCutsCouples();

        auto material = G4Material::GetMaterial("GAGG");
        if (!material) {
            throw std::runtime_error("Cannot find material GAGG.");
        }

        const auto energyGrid = MakeLogEnergyGrid(
            Settings::kEminMeV,
            Settings::kEmaxMeV,
            Settings::kNPoints
        );

        std::cout << "Generating stopping-power tables for GAGG" << std::endl;
        std::cout << "  E range : " << Settings::kEminMeV
                  << " to " << Settings::kEmaxMeV << " MeV" << std::endl;
        std::cout << "  points  : " << Settings::kNPoints << std::endl;
        std::cout << "  density : " << material->GetDensity()/(g/cm3)
                  << " g/cm3" << std::endl;

        WriteStoppingTable(
            Settings::kElectronOut,
            "e-",
            "eIoni",
            material,
            energyGrid
        );

        WriteStoppingTable(
            Settings::kAlphaOut,
            "alpha",
            "ionIoni",
            material,
            energyGrid
        );

        // 16O recoil/fragment stopping table.
        // This table uses the total kinetic energy of the O16 ion, not MeV/u.
        // It can be used to build an O16-specific Birks yield curve if
        // alphaInelastic creates O16 tracks in GAGG.
        auto oxygen16 = GetIonDefinition(8, 16, 0.0);
        WriteStoppingTableForParticle(
            Settings::kO16Out,
            oxygen16,
            "O16",
            "ionIoni",
            material,
            energyGrid
        );
    }
};


int main()
{
    auto runManager = new G4RunManager();

    runManager->SetUserInitialization(new DetectorConstruction());

    G4PhysListFactory factory;
    auto physicsList = factory.GetReferencePhysList("FTFP_BERT_EMZ");
    physicsList->SetDefaultCutValue(1.0*um);
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserAction(new RunAction());

    runManager->Initialize();

    // BeamOn(1) is used only to enter the normal run lifecycle.
    // The CSV files are written in BeginOfRunAction, before the dummy event.
    runManager->BeamOn(1);

    delete runManager;
    return 0;
}
