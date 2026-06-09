# MCGPU_materials

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andreubs/MCGPU_materials/HEAD?labpath=create_MCGPU_material.ipynb)

Simple utility to generate `.mcgpu` material files for the **MC-GPU** GPU-accelerated 
Monte Carlo X-ray transport code without a local PENELOPE 2006 installation 
(material files compatible with versions v1.3, VICTRE_MCGPU, MCGPU-PET). 
A single interactive Fortran program replaces the original two-step
workflow (`material.f` -> `.mat` -> `MC-GPU_create_material_data.f` -> `.mcgpu`).

---

## Quickstart: run in the browser

Click the **Binder** badge above. The notebook `create_MCGPU_material.ipynb`
opens with gfortran pre-installed and the binary pre-compiled. Edit the material
parameters in Step 2, run all cells, and download your `.mcgpu` file.
No local installation required.

---

## Local use

### Prerequisites

- gfortran (GCC Fortran, version 5 or later)
- Python 3 with numpy and matplotlib (only for the notebook and validation script)
- The `PENDBASE_photons/` database is included in this repository -- no setup needed.

### Compile

```
gfortran MCGPU_materials.f penelope_photons.f -o MCGPU_materials.x -O3
```

### Run

```
./MCGPU_materials.x
```

or non-interactively with a redirect file:

```
./MCGPU_materials.x < my_material.txt
```

### Interactive prompts (in order)

The program asks the user for the information required to define the materials.
These questions where originally requested by PENELOPE's `material.f` and
MC-GPU's `MC-GPU_create_material_data.f` (see their documentation for reference):

1. **Mode**: `1` keyboard entry, `2` look up by ID from `PENDBASE_photons/pdcompos.p06`
   (IDs 1-99 = elements, 100-280 = compounds)

2. **Composition** (keyboard mode): material name, number of elements, then for
   each element Z and stoichiometric index (atoms/molecule); then mass density (g/cm3).
   Single-element materials skip the stoichiometric format prompt.

3. **Insulator?** `1` = yes (polymers, tissue, glass -- sets Fcb=Wcb=0).
   `2` = no/auto (metals and semiconductors -- free-electron plasmon model for
   outer-shell electrons). Affects the Compton shell structure; negligible effect
   on total MFP at diagnostic energies.

4. **Energy range**: `Emin Emax` in eV. One extra bin above Emax is added
   automatically so MC-GPU can interpolate up to exactly Emax.

5. **Energy step DE** in eV. Choose DE so that `(Emax-Emin)/DE` is an integer
   to get exact round-number energies (e.g. DE=5 for 5-120 keV).
   Number of bins = `NINT((Emax-Emin)/DE) + 2`.

6. **Output filename** (e.g. `water_5-120keV.mcgpu`)

### Example input file (liquid water, 5-120 keV)

The following text file can be redirected to the executable stdin instead of manually typing the information:

```
1
Water, liquid
2
1 2
8 1
1.0
1
5000 120000
5
water_5-120keV.mcgpu
```

---

## Quick reference: common materials

### By keyboard entry

| Material | ELEMENTS (Z stoich) | Density (g/cm3) | Insulator |
|----------|---------------------|-----------------|-----------|
| Water H2O | `1 2` then `8 1` | 1.000 | yes |
| PMMA C5H8O2 | `6 5`, `1 8`, `8 2` | 1.19 | yes |
| Polyethylene (CH2)n | `6 1` then `1 2` | 0.94 | yes |
| Aluminum | `13` | 2.699 | no |
| Copper | `29` | 8.96 | no |
| Tungsten | `74` | 19.3 | no |
| CsI (scintillator) | `55 1` then `53 1` | 4.51 | yes |

### From the PENELOPE database (mode 2)

| ID | Material |
|----|----------|
| 104 | Air, dry |
| 119 | Bone, cortical (ICRU) |
| 153 | Fat, adipose tissue |
| 155 | Glandular tissue |
| 174 | Lung tissue |
| 194 | Muscle, skeletal |
| 216 | PMMA |
| 226 | Polyethylene |
| 243 | Silicon |
| 276 | Polycarbonate |
| 278 | Water, liquid |

Full list: `PENDBASE_photons/pdcompos.p06` (material name on the first line of
each entry).

---

## Output file format (.mcgpu)

Three sections written sequentially:

**1 -- MFP table** (one row per linear energy bin)

```
Energy(eV)   MFP_Rayleigh   MFP_Compton   MFP_Photoelectric   MFP_Total   Rayleigh_Pmax
```

MFP_Total includes pair production for correct beam attenuation above 1.022 MeV.
Below 1.022 MeV pair production contributes zero.

**2 -- Rayleigh RITA sampling block** (128 rows)

```
X   P   A   B   ITL   ITU
```

Adaptive grid for sampling coherent scattering angles from the squared molecular
form factor F2(x). Analytical Balyuzi form factors from PENELOPE BLOCK DATA PENDAT.

**3 -- Compton shell block** (one row per electron shell group)

```
FCO   UICO   FJ0   KZCO   KSCO
```

| Column | Meaning |
|--------|---------|
| FCO | Electrons in shell group |
| UICO | Shell ionisation energy (eV) |
| FJ0 | One-electron Compton profile at pz=0 (scaled by 1/alpha) |
| KZCO | Atomic number of element owning this shell (0 = conduction band) |
| KSCO | Shell index (1=K, 2=L1, ..., 30=outer/grouped) |

Format: `1pe13.6` (6 significant figures). Compatible with MC-GPU v1.3,
VICTRE_MCGPU, and MCGPU-PET.

---

## Database: PENDBASE_photons

200-file photon-only subset of the PENELOPE 2006 PENDBASE (797 files, ~160 MB).
Included in this repository -- no separate download needed.

| Files | Count | Content |
|-------|-------|---------|
| `pdgph01.p06` ... `pdgph99.p06` | 99 | Photoelectric cross sections (Z=1-99) |
| `pdgpp01.p06` ... `pdgpp99.p06` | 99 | Pair production cross sections (Z=1-99) |
| `pdatconf.p06` | 1 | Atomic shell structure and Compton profiles |
| `pdcompos.p06` | 1 | 280 pre-defined material compositions |

Electron/positron interaction files (~140 MB) are not included -- not needed for
photon-only transport.

---

## Physics

Cross sections come from the PENELOPE 2006 library via `penelope_photons.f`, a
1377-line photon-only subset (16 subroutines) of the full `penelope.f` library.

| Interaction | Source |
|-------------|--------|
| Rayleigh | Analytical Balyuzi form factors (hardcoded in BLOCK DATA PENDAT) |
| Compton | Relativistic impulse approximation with one-electron profiles from pdatconf.p06 |
| Photoelectric | Tabulated cross sections from pdgph##.p06 |
| Pair production | Tabulated cross sections from pdgpp##.p06. For energies above MeV, the pair production attenuation is added to the total mean free path, but no e+/e- are generated. |

Known limitations:
- Photoelectrons and fluorescence not included (MC-GPU does not track secondary particles)
- Pair production attenuates the beam correctly but generates no secondary particles
- Single material per file

---

## Validation

The `material_validation/` folder contains reference `.mcgpu` files for aluminum
generated by the original PENELOPE 2006 two-step workflow, and a Python script
`test_aluminum_compatibility.py` that compares them against files produced by
the new tool.

```
cd material_validation
python3 test_aluminum_compatibility.py
```

Expected results (aluminum, two energy ranges):

| Section | Max relative difference |
|---------|------------------------|
| Compton MFP | < 0.02% |
| Photoelectric MFP | < 0.001% |
| Total MFP | < 0.004% |
| Rayleigh MFP | < 0.11% (known; direct evaluation vs. PENELOPE spline) |
| Rayleigh RITA block | < 0.001% |
| Compton shells (FJ0) | < 0.001% |
| Integer fields (ITL, ITU, KZCO, KSCO) | exact match |

---

## Binder configuration

The `binder/` folder contains files required for the Binder launch badge.
The files request install of gfortran and other dependencies and compile the code.

---

## References

1. F. Salvat, J. M. Fernandez-Varea and J. Sempau,
   *PENELOPE -- A code system for Monte Carlo simulation of electron and photon
   transport*, NEA-OECD, Issy-les-Moulineaux (2006)

2. A. Badal and A. Badano,
   *Accelerating Monte Carlo simulations of photon transport in a voxelized
   geometry using a massively parallel Graphics Processing Unit*,
   Med. Phys. **36**, pp. 4878-4880 (2009)

3. MC-GPU repositories:
   - v1.3: https://github.com/DIDSR/MCGPU
   - VICTRE_MCGPU: https://github.com/DIDSR/VICTRE_MCGPU
   - MCGPU-PET: https://github.com/DIDSR/MCGPU-PET

---

## PENELOPE copyright notice (penelope.f and PENDBASE)

```                                                         
  PENELOPE/PENGEOM (version 2006)                                     
  Copyright (c) 2001-2006                                             
  Universitat de Barcelona                                            

  Permission to use, copy, modify, distribute and sell this software  
  and its documentation for any purpose is hereby granted without     
  fee, provided that the above copyright notice appears in all        
  copies and that both that copyright notice and this permission      
  notice appear in all supporting documentation. The Universitat de   
  Barcelona makes no representations about the suitability of this    
  software for any purpose. It is provided "as is" without express    
  or implied warranty.                                                
```

---
