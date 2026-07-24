# GAGGWaveformSim

Minimal Geant4 11.3.2 project for GAGG optical-photon waveform truth output.

Geometry in the first version:

- GAGG crystal: 10 mm x 10 mm x 1 mm
- Incident face: 10 mm x 10 mm, source enters from the -z side
- Virtual SiPM layer: 10 mm x 10 mm x 0.01 mm, attached to the +z side

The virtual SiPM layer is a sensitive absorbing optical-photon plane. It records photon arrival times and kills optical photons after recording. It does not model PDE, avalanche, cross talk, afterpulse, or electronics.

ROOT output:

- RunInfo
- WaveformEvent
- PrimaryPhoton

Time unit for photon arrival lists is ps. `UInt_t` is used for photon time, so the current 1 ms window corresponds to 0 to 1,000,000,000 ps.

Build:

```bash
mkdir build
cd build
cmake ..
make -j
./GAGGWaveformSim ../macros/run.mac
```

Output ROOT files are written to `./data` under the build directory.
