# material_validation

Reference .mcgpu files and compatibility test for aluminum.

## Reference files

| File | Material | Energy range | Step | Bins | Origin |
|------|----------|-------------|------|------|--------|
| `Al_120keV_REFERENCE.mcgpu` | Aluminum (Z=13, 2.699 g/cm3) | 5-120 keV | 5 eV | 23002 | Legacy PENELOPE 2006 two-step workflow |
| `Al_10MeV_REFERENCE.mcgpu` | Aluminum (Z=13, 2.699 g/cm3) | 100 keV-10 MeV | 100 keV | 100 | Legacy PENELOPE 2006 two-step workflow |

The reference files were produced by the original two-step workflow:
`material.f` (PENELOPE 2006) -> `.mat` file -> `MC-GPU_create_material_data.f`.
They use 10-significant-figure output (`1pe17.10` format).

## Validation test

`test_aluminum_compatibility.py` compares the GENERATED files against the REFERENCE
files and reports the maximum relative difference for every column in all three
sections of the .mcgpu format (MFP table, Rayleigh RITA block, Compton shells).

The GENERATED files must be produced first:

```bash
gfortran ../MCGPU_materials.f ../penelope_photons.f -o ../MCGPU_materials.x -O3
../MCGPU_materials.x < Al_120keV_input.txt
../MCGPU_materials.x < Al_10MeV_input.txt
python3 test_aluminum_compatibility.py
```

Expected results:
- All columns: max relative difference < 0.05%
- Rayleigh MFP (diagnostic case only): max ~0.10% -- known difference between
  direct cross-section evaluation (new tool) and PENELOPE spline interpolation
  (legacy tool). Both are physically correct; the threshold is set to 0.20%.

The script also saves `aluminum_validation_plots.png` with a 2x2 figure:
MFP curves (reference solid, generated open circles) and relative differences.

## Input files

| File | Description |
|------|-------------|
| `Al_120keV_input.txt` | Input for the 5-120 keV aluminum case |
| `Al_10MeV_input.txt` | Input for the 100 keV-10 MeV aluminum case |
| `polycarbonate_input.txt` | Example input for polycarbonate C15H16O2 |
