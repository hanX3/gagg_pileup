#ifndef Constants_H
#define Constants_H 1

#include <cstdint>
#include <cmath>

#include "G4SystemOfUnits.hh"

// ROOT output directory. Run from build/; data are written outside build/.
static const char DATAPATH[] = "../data";

// One Geant4 event is one waveform time window.
// For pileup studies this is both the gamma-emission window and the recorded
// photon-arrival window.
static const std::uint32_t EVENT_TIME_LENGTH_PS = 10000000U; // 10 us

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
// Primary source defaults
// -----------------------------------------------------------------------------
// The simple source used for alpha/gamma waveform comparison is a circular
// surface source centered at SOURCE_POSITION_* and emitting into a cone around
// +z toward the GAGG.  The cone is intentionally narrow to avoid wasting events
// on particles that cannot geometrically reach the detector region.
static const bool DEFAULT_USE_DISK_CONE_SOURCE = false;

// phi20 mm surface source.
static const double SOURCE_DISK_DIAMETER_MM = 20.0;
static const double SOURCE_DISK_RADIUS_MM = 0.5*SOURCE_DISK_DIAMETER_MM;

// Cone-angle estimate.  The target radius is taken as source radius plus the
// GAGG half diagonal, then enlarged by a safety factor.
static const double PI_VALUE = 3.1415926535897932384626433832795;
static const double GAGG_HALF_DIAGONAL_MM =
  0.5*std::sqrt(GAGG_SIZE_X_MM*GAGG_SIZE_X_MM + GAGG_SIZE_Y_MM*GAGG_SIZE_Y_MM);
static const double SOURCE_TO_GAGG_FRONT_DISTANCE_MM =
  GAGG_FRONT_SURFACE_Z_MM - SOURCE_POSITION_Z_MM;
static const double SOURCE_EMISSION_CONE_SAFETY_FACTOR = 1.20;
static const double SOURCE_EMISSION_CONE_TARGET_RADIUS_MM =
  SOURCE_EMISSION_CONE_SAFETY_FACTOR*(SOURCE_DISK_RADIUS_MM + GAGG_HALF_DIAGONAL_MM);
static const double SOURCE_EMISSION_CONE_HALF_ANGLE_DEG =
  std::atan(SOURCE_EMISSION_CONE_TARGET_RADIUS_MM/SOURCE_TO_GAGG_FRONT_DISTANCE_MM)*180.0/PI_VALUE;

// Aim-at-GAGG source mode.  In this mode, each event samples a random point on
// the phi20-mm source disk and a random target point on the GAGG front face.
// The primary momentum direction is the line from source point to target point.
// This is an efficient conditional source for waveform studies: it samples
// events that enter the GAGG, not the absolute geometrical efficiency of a real
// isotropic source.
static const bool DEFAULT_AIM_AT_GAGG = false;
static const double SOURCE_TARGET_CENTER_X_MM = 0.0;
static const double SOURCE_TARGET_CENTER_Y_MM = 0.0;
static const double SOURCE_TARGET_Z_MM = GAGG_FRONT_SURFACE_Z_MM;
static const double SOURCE_TARGET_SIZE_X_MM = GAGG_SIZE_X_MM;
static const double SOURCE_TARGET_SIZE_Y_MM = GAGG_SIZE_Y_MM;


// -----------------------------------------------------------------------------
// Pileup particle/photon-source defaults
// -----------------------------------------------------------------------------
// One Geant4 event corresponds to one waveform time window.  In pileup mode,
// each enabled source component is sampled independently with Poisson(rate*T),
// where T = EVENT_TIME_LENGTH_PS in seconds.  All generated primaries are
// then sorted by emission time, so primary_id follows chronological order.
static const bool DEFAULT_PILEUP_ENABLED = false;

// p-11B-related source components used in the current parameterized model:
//   1. bremsstrahlung: p-11B plasma X-ray/photon background;
//   2. c12_capture: 11B(p,gamma)12C capture/de-excitation photon component;
//   3. alpha_p11b: effective alpha-particle source from p + 11B -> 3 alpha.
//
// The 11B(p,p'gamma)11B 2.125-MeV inelastic-gamma component has been removed
// from this model because it is not expected in the energy range currently
// considered for the p-11B fusion case.
//
// Only component enable flags and rates are intended to be changed in macros.
static const bool DEFAULT_BREMS_SOURCE_ENABLED = false;
static const bool DEFAULT_C12_CAPTURE_SOURCE_ENABLED = false;
static const bool DEFAULT_ALPHA_P11B_SOURCE_ENABLED = false;

static const double DEFAULT_BREMS_SOURCE_RATE_HZ = 1.0e4;
static const double DEFAULT_C12_CAPTURE_SOURCE_RATE_HZ = 1.0e3;
static const double DEFAULT_ALPHA_P11B_SOURCE_RATE_HZ = 1.0e4;

// Maxwellian thermal bremsstrahlung photon-number spectrum:
//   dN/dE proportional to E^{-1} exp(-E/kTe), Emin <= E <= Emax.
// These values define the nominal p-11B plasma photon background.  The nominal
// kTe is 100 keV; 40--300 keV can be used in dedicated code-level scans if
// needed, but it is intentionally not exposed as a routine macro parameter here.
static const double BREMS_EMIN = 10.0*keV;
static const double BREMS_EMAX = 5.0*MeV;
static const double BREMS_KTE  = 100.0*keV;


// 11B(p,gamma)12C capture/de-excitation photon component.
// Reference for low-energy 11B(p,gamma)12C gamma branches:
//   J. J. He et al., "Direct measurement of 11B(p,gamma)12C astrophysical
//   S factors at low energies", Phys. Rev. C 93, 055804 (2016).
// In that notation:
//   gamma0      : capture to 12C ground state, E_gamma ~ 16.1 MeV;
//   gamma1      : capture to 12C first excited state, E_gamma ~ 11.7 MeV;
//   gamma_star  : 12C*(4.439) -> 12C(g.s.) de-excitation gamma.
// The measured low-energy yield ratio gamma0/gamma1 is approximately 4.6%.
static const double C12_CAPTURE_GAMMA0_ENERGY = 16.1*MeV;
static const double C12_CAPTURE_GAMMA1_ENERGY = 11.7*MeV;
static const double C12_CAPTURE_GAMMA_STAR_ENERGY = 4.439*MeV;
static const double C12_CAPTURE_GAMMA0_TO_GAMMA1_RATIO = 0.046;

// Convert gamma0/gamma1 yield ratio to reaction-branch probabilities:
//   B0 = gamma0 branch, B1 = gamma1 branch, B0/B1 = 0.046.
static const double C12_CAPTURE_BRANCH_GROUND =
  C12_CAPTURE_GAMMA0_TO_GAMMA1_RATIO /
  (1.0 + C12_CAPTURE_GAMMA0_TO_GAMMA1_RATIO);
static const double C12_CAPTURE_BRANCH_FIRST_EXCITED =
  1.0 / (1.0 + C12_CAPTURE_GAMMA0_TO_GAMMA1_RATIO);

// Current implementation uses an effective single-photon line-mixture model for
// the c12_capture component.  It does not explicitly force gamma1 and
// gamma_star to be generated as a coincidence pair.  The photon-yield weights
// below follow from the reaction branches:
//   gamma0 branch emits one photon;
//   gamma1 branch emits gamma1 plus gamma_star, i.e. two photons.
static const double C12_CAPTURE_PHOTON_WEIGHT_DENOMINATOR =
  C12_CAPTURE_BRANCH_GROUND + 2.0*C12_CAPTURE_BRANCH_FIRST_EXCITED;
static const double C12_CAPTURE_GAMMA0_PHOTON_WEIGHT =
  C12_CAPTURE_BRANCH_GROUND / C12_CAPTURE_PHOTON_WEIGHT_DENOMINATOR;
static const double C12_CAPTURE_GAMMA1_PHOTON_WEIGHT =
  C12_CAPTURE_BRANCH_FIRST_EXCITED / C12_CAPTURE_PHOTON_WEIGHT_DENOMINATOR;
static const double C12_CAPTURE_GAMMA_STAR_PHOTON_WEIGHT =
  C12_CAPTURE_BRANCH_FIRST_EXCITED / C12_CAPTURE_PHOTON_WEIGHT_DENOMINATOR;


// p + 11B -> 3 alpha effective single-alpha source spectrum.
// Current model is used for detector-response and pileup studies, not for a
// full three-body Dalitz simulation.  The source rate is interpreted as an
// effective alpha-particle rate toward the detector, not as a fusion-reaction
// rate.
//
// Model:
//   1/3 primary alpha:   E = 3.76 MeV;
//   2/3 secondary alpha: sequential 8Be* decay approximation,
//      E = Eboost + Estar + 2*sqrt(Eboost*Estar)*cos(theta),
//      cos(theta) uniformly sampled in [-1, 1].
static const double P11B_ALPHA_PRIMARY_FRACTION = 1.0/3.0;
static const double P11B_ALPHA_PRIMARY_ENERGY = 3.76*MeV;
static const double P11B_ALPHA_SECONDARY_ESTAR = 1.515*MeV;
static const double P11B_ALPHA_SECONDARY_EBOOST = 0.94*MeV;

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

// Two-component scintillation time constants used by Geant4 default
// G4Scintillation.  The current model uses the same decay constants for
// gamma/electron and alpha responses; particle dependence is introduced only
// through the total yield curves and the component fractions below.
static const double GAGG_FAST_TIME_NS = 95.0;
static const double GAGG_SLOW_TIME_NS = 351.0;

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

// Alpha pulse-shape component fractions.  These are used directly in
// DetectorConstruction.cc as Geant4 particle-dependent component fractions.
static const double ALPHA_FAST_FRACTION = 0.30;
static const double ALPHA_SLOW_FRACTION = 0.70;

// Generic heavy-ion / recoil-ion response parameters.
// Geant4 uses the generic IONSCINTILLATIONYIELD key for heavy ions and
// recoil nuclei such as O16, Al, Ga, and Gd.  The corresponding total
// light-yield vector should be named gaggIonScintillationYield in
// gagg_birks_yield_curves.inc.  In the current model, its two-component
// pulse-shape fraction is treated as alpha-like because no independent
// heavy-ion PSD fraction is used.
static const double ION_FAST_FRACTION = ALPHA_FAST_FRACTION;
static const double ION_SLOW_FRACTION = ALPHA_SLOW_FRACTION;

// Backward-compatible name for total gamma/electron scintillation yield.
// The current Geant4 material table uses the externally generated cumulative
// yield curves instead of a single global SCINTILLATIONYIELD constant.
static const double GAGG_SCINT_YIELD_PER_MEV = GAMMA_ELECTRON_LIGHT_YIELD_PER_MEV;

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
