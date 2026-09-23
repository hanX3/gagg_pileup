#ifndef Constants_H
#define Constants_H 1

#include <cstdint>

// ROOT output directory. Run from build/; data are written outside build/.
static const char DATAPATH[] = "../data";

// One Geant4 event is one waveform time window.
static const std::uint32_t EVENT_TIME_LENGTH_PS = 1000000000U; // 1 ms

// -----------------------------------------------------------------------------
// Geometry constants, in mm unless otherwise stated
// -----------------------------------------------------------------------------
// World is a 1 m cube centered at the origin.
static const double WORLD_SIZE_MM = 1000.0;

// GAGG crystal
static const double GAGG_DENSITY_G_CM3 = 6.63;
static const double GAGG_SIZE_X_MM = 10.0;
static const double GAGG_SIZE_Y_MM = 10.0;
static const double GAGG_SIZE_Z_MM = 1.0;

// Virtual SiPM sensitive layer
static const double SIPM_SIZE_X_MM = GAGG_SIZE_X_MM;
static const double SIPM_SIZE_Y_MM = GAGG_SIZE_Y_MM;
static const double SIPM_SIZE_Z_MM = 0.01;

// Front entrance Mylar window: 5 um = 0.005 mm
static const double MYLAR_FRONT_SIZE_X_MM = GAGG_SIZE_X_MM;
static const double MYLAR_FRONT_SIZE_Y_MM = GAGG_SIZE_Y_MM;
static const double MYLAR_FRONT_SIZE_Z_MM = 0.005;

// Four side Mylar films: 20 um = 0.020 mm
static const double MYLAR_SIDE_THICKNESS_MM = 0.020;

// Particle source position.
// Current simple particle gun emits from the origin and points along +z.
static const double SOURCE_POSITION_X_MM = 0.0;
static const double SOURCE_POSITION_Y_MM = 0.0;
static const double SOURCE_POSITION_Z_MM = 0.0;

// Detector position along +z.
// Incoming particles travel from the origin along +z.
// GAGG entrance surface = -z surface of GAGG, placed at z = 250 mm.
static const double GAGG_FRONT_SURFACE_Z_MM = 250.0;
static const double GAGG_CENTER_Z_MM = GAGG_FRONT_SURFACE_Z_MM + 0.5*GAGG_SIZE_Z_MM;
static const double SIPM_CENTER_Z_MM = GAGG_FRONT_SURFACE_Z_MM + GAGG_SIZE_Z_MM + 0.5*SIPM_SIZE_Z_MM;
static const double MYLAR_FRONT_CENTER_Z_MM = GAGG_FRONT_SURFACE_Z_MM - 0.5*MYLAR_FRONT_SIZE_Z_MM;

// -----------------------------------------------------------------------------
// Optical constants
// -----------------------------------------------------------------------------
static const double VACUUM_RINDEX = 1.00;
static const double GAGG_RINDEX = 1.90;
static const double VIRTUAL_SIPM_RINDEX = GAGG_RINDEX;
static const double MYLAR_RINDEX = 1.65;

static const double GAGG_ABS_LENGTH_CM = 100.0;
static const double MYLAR_ABS_LENGTH_CM = 100.0;

// -----------------------------------------------------------------------------
// GAGG scintillation and response-label parameters
// -----------------------------------------------------------------------------
// Current analysis scheme:
//   edep_total_MeV is the Geant4 truth energy deposited in GAGG.
//   optical-photon counts and arrival times are the detector optical signal.
//   No independent light_output_MeVee branch is stored.
static const double GAMMA_ELECTRON_LIGHT_YIELD_PER_MEV = 46000.0;
static const double GAMMA_ELECTRON_FAST_TIME_NS = 95.0;
static const double GAMMA_ELECTRON_SLOW_TIME_NS = 351.0;
static const double GAMMA_ELECTRON_FAST_FRACTION = 0.79;
static const double GAMMA_ELECTRON_SLOW_FRACTION = 0.21;

// Geant4 Birks constant derived from the GAGG(Ce) LET coefficient.
// It affects Geant4 scintillation-photon generation directly.
// LET_mass coefficient unit: g cm-2 MeV-1.
static const double ALPHA_BIRKS_A1_G_CM2_PER_MEV = 6.5e-3;

// Equivalent Geant4 Birks constant kB in mm/MeV, derived from
// ALPHA_BIRKS_A1 / density.
static const double GAGG_BIRKS_CONSTANT_MM_PER_MEV =
  (ALPHA_BIRKS_A1_G_CM2_PER_MEV/GAGG_DENSITY_G_CM3)*10.0;

// First alpha pulse-template parameters for later waveform/post-processing.
// The Geant4 material scintillation table below still uses the global
// gamma/electron time constants unless a custom scintillation model is added.
static const double ALPHA_FAST_TIME_NS = 95.0;
static const double ALPHA_SLOW_TIME_NS = 351.0;
static const double ALPHA_FAST_FRACTION = 0.30;
static const double ALPHA_SLOW_FRACTION = 0.70;

// Backward-compatible names used by DetectorConstruction.cc.
// They intentionally correspond to the gamma/electron global scintillation model.
static const double GAGG_SCINT_YIELD_PER_MEV = GAMMA_ELECTRON_LIGHT_YIELD_PER_MEV;
static const double GAGG_FAST_TIME_NS = GAMMA_ELECTRON_FAST_TIME_NS;
static const double GAGG_SLOW_TIME_NS = GAMMA_ELECTRON_SLOW_TIME_NS;
static const double GAGG_FAST_FRACTION = GAMMA_ELECTRON_FAST_FRACTION;
static const double GAGG_SLOW_FRACTION = GAMMA_ELECTRON_SLOW_FRACTION;

// -----------------------------------------------------------------------------
// Step-limit defaults for charged-particle transport in sensitive volumes
// -----------------------------------------------------------------------------
// These values are used by G4UserLimits. With G4StepLimiterPhysics default
// behavior, the limit applies to charged particles only, not to gamma rays or
// optical photons. They can be overridden in the macro before /run/initialize.
static const double DEFAULT_GAGG_MAX_STEP_UM = 1.0;
static const double DEFAULT_MYLAR_FRONT_MAX_STEP_UM = 0.5;
static const double DEFAULT_MYLAR_SIDE_MAX_STEP_UM = 1.0;

#endif
