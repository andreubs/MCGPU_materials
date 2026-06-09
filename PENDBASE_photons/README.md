# PENDBASE_photons

Photon-only subset of the PENELOPE 2006 PENDBASE cross-section database, containing the 200 files needed by `MCGPU_material_generator` to generate `.mcgpu` material files for MC-GPU photon transport simulations.

## Contents

| Files | Count | Description |
|-------|-------|-------------|
| `pdgph01.p06` … `pdgph99.p06` | 99 | Photoelectric effect cross sections (one file per element, Z = 1–99) |
| `pdgpp01.p06` … `pdgpp99.p06` | 99 | Pair production cross sections (one file per element, Z = 1–99) |
| `pdatconf.p06` | 1 | Atomic electron shell structure, binding energies, and one-electron Compton profiles |
| `pdcompos.p06` | 1 | 280 pre-defined material compositions (elements and compounds of radiological interest) |

> **Note on pair production:** MC-GPU does not simulate pair production events, but these cross sections are included in the total attenuation coefficient for accurate photon beam fluence at energies above the pair production threshold (1.022 MeV). At diagnostic imaging energies (< 150 keV) their contribution is exactly zero.

## What was removed from the full PENDBASE (797 files)

The ~140 MB of electron/positron interaction files are not needed for photon-only transport:

| Files | Content |
|-------|---------|
| `pdeel##.p06` | Electron elastic scattering |
| `pdesi##.p06` | Electron inner-shell ionization by electron impact |
| `pdpsi##.p06` | Inner-shell ionization by positron impact |
| `pdebr##.p06` | Electron bremsstrahlung emission data |
| `peldx###.p06` | Electron elastic differential cross sections (~700 KB/element) |
| `eeldx###.p06` | Positron elastic differential cross sections (~700 KB/element) |
| `pdrelax.p06` | Atomic relaxation cascade data (fluorescence X-rays) — MC-GPU does not track secondary particles |
| `pdbrang.p06` | Bremsstrahlung angular distribution (electron-specific) |
| `pdeflist.p06` | Material definition list (informational only) |

## Reference

- F. Salvat, J. M. Fernandez-Varea and J. Sempau, *PENELOPE — A code system for Monte Carlo simulation of electron and photon transport*, NEA-OECD, Issy-les-Moulineaux (2006)
