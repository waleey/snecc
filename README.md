# SNECC: Supernova Neutrino Energy Constraints Calculator

This repository contains a command-line pipeline for running the IceCube supernova neutrino coincidence-ratio analysis.

The code takes supernova neutrino flux models, applies an oscillation scenario, folds the model through detector-hit information and generated spectra, and saves the coincidence-ratio results and combined chi-square confidence contours for one or more supernova distances.

For each requested model, mass hierarchy, and distance, the script saves:

- coincidence-ratio output data,
- combined chi-square contour data,
- 68% and 99% confidence intervals,
- a JSON input deck recording the run configuration.

---

## Clone the Repository

Clone this repository with:

```bash
git clone **https://github.com/waleey/snecc.git**
cd **snecc**
```

---

## Main Python Script

The main analysis script is:

```text
snecc_main_optimized.py
```

This script is normally run through the provided Bash script:

```text
run_snecc.sh
```

---

## Input Variables

The main Python script accepts the following command-line arguments.

### `--modelName`

Name of the model being analyzed.

This name is used in the output directory structure.

Example:

```bash
--modelName sukhbold_sfho
```

Users should change this if running a different model set.

Change this:

**`sukhbold_sfho`**

---

### `--modelPath`

Path to the SNEWPY-compatible model file or directory.

This can be either:

1. A single `.h5` model file, or
2. A directory containing multiple `.h5` model files.

Example:

```bash
--modelPath /Users/walu/icecube/snewpy_model_data/sukhbold_sfho/
```

Users must change this path to point to their own model file or model directory.

Change this:

**`/Users/walu/icecube/snewpy_model_data/sukhbold_sfho/`**

---

### `--minMultiplicity`

Minimum coincidence multiplicity included in the analysis.

Example:

```bash
--minMultiplicity 2
```

---

### `--maxMultiplicity`

Maximum coincidence multiplicity included in the analysis.

Example:

```bash
--maxMultiplicity 15
```

---

### `--massHierarchy`

Neutrino oscillation scenario.

Allowed options are:

```text
nmo
imo
no_osc
all
```

where:

- `nmo` means normal mass ordering,
- `imo` means inverted mass ordering,
- `no_osc` means no oscillation,
- `all` runs all three options.

Example:

```bash
--massHierarchy nmo
```

Change this if you want a different oscillation scenario:

**`nmo`**

---

### `--nGen`

Number of generated neutrinos/events used in the detector simulation using Geant4.

Example:

```bash
--nGen 5255019
```

Change this to match your generated simulation sample:

**`5255019`**

---

### `--distance`

One or more supernova distances in kiloparsecs.

Example:

```bash
--distance 3.0 5.0 10.0 15.0 20.0 25.0
```

The script loops over every distance provided here.

Change these distances as needed:

**`3.0 5.0 10.0 15.0 20.0 25.0`**

---

### `--nModules`

Number of detector modules used in the analysis.

Example:

```bash
--nModules 15000
```

Change this if your detector configuration uses a different number of modules:

**`15000`**

---

### `--hitFilePath`

Path to the detector-hit file.

Example:

```bash
--hitFilePath /Users/walu/icecube/doumeki_analysis/mdom_ibd_combined_52550199_nu_energy_flat_spectrum.csv
```

This file should contain the simulated detector-hit information generated using Geant4 simulations.

Change this path:

**`/Users/walu/icecube/doumeki_analysis/mdom_ibd_combined_52550199_nu_energy_flat_spectrum.csv`**

---

### `--spectraFilePath`

Path to the spectra grid file.

Example:

```bash
--spectraFilePath /Users/walu/icecube/energy_constraints/snecc/data/spectra_grid.h5
```

Change this path:

**`/Users/walu/icecube/energy_constraints/snecc/data/spectra_grid.h5`**

---

### `--output_dir`

Base directory where output files will be saved.

Example:

```bash
--output_dir output/
```

Change this if you want the results saved somewhere else:

**`output/`**

---

### `--bkgFilePath`

Path to the background multiplicity file.

This argument is optional. If it is not provided, the script uses the default path:

```bash
../../OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv
```

To provide your own background file, add:

```bash
--bkgFilePath **/path/to/background_file.csv**
```

Change this if your background file is somewhere else:

**`/path/to/background_file.csv`**

---

## How to Run

The recommended way to run the analysis is through the provided Bash script.

First, make the script executable:

```bash
chmod +x run_snecc.sh
```

Then run:

```bash
./run_snecc.sh
```

---

## Provided Bash Script

The included Bash script runs the main Python script from the command line:

```bash
#!/bin/bash

# ============================================================
# Supernova Energy Constraint Big Run Script
# ============================================================

set -e
set -o pipefail

echo "============================================================"
echo "Starting SNECC run..."
echo "============================================================"

CMD="python3 snecc_main_optimized.py \
  --modelName sukhbold_sfho \
  --modelPath /Users/walu/icecube/snewpy_model_data/sukhbold_sfho/ \
  --minMultiplicity 2 \
  --maxMultiplicity 15 \
  --massHierarchy nmo \
  --nGen 5255019 \
  --distance 3.0 5.0 10.0 15.0 20.0 25.0 \
  --nModules 15000 \
  --hitFilePath /Users/walu/icecube/doumeki_analysis/mdom_ibd_combined_52550199_nu_energy_flat_spectrum.csv \
  --spectraFilePath /Users/walu/icecube/energy_constraints/snecc/data/spectra_grid.h5 \
  --output_dir output/"

echo "Executing:"
echo "$CMD"
echo "------------------------------------------------------------"

eval $CMD

echo "============================================================"
echo "Run completed successfully."
echo "============================================================"
```

Before running it, update the machine-specific paths:

- **`/Users/walu/icecube/snewpy_model_data/sukhbold_sfho/`**
- **`/Users/walu/icecube/doumeki_analysis/mdom_ibd_combined_52550199_nu_energy_flat_spectrum.csv`**
- **`/Users/walu/icecube/energy_constraints/snecc/data/spectra_grid.h5`**
- **`output/`**

Also update these values if needed:

- **`sukhbold_sfho`**
- **`nmo`**
- **`5255019`**
- **`3.0 5.0 10.0 15.0 20.0 25.0`**
- **`15000`**
- **`2`**
- **`15`**

---

## Output Structure

The output directory is organized as:

```text
output_dir/
└── modelName/
    └── massHierarchy/
        └── distance/
            ├── *_ratio.h5
            ├── *_chisq.h5
            └── *_inputdeck.json
```

For the provided Bash script, the output will look like:

```text
output/
└── sukhbold_sfho/
    └── nmo/
        ├── 3kpc/
        ├── 5kpc/
        ├── 10kpc/
        ├── 15kpc/
        ├── 20kpc/
        └── 25kpc/
```

Inside each distance folder, the code saves files with names like:

```text
MODELFILE_nmo_10kpc_ratio.h5
MODELFILE_nmo_10kpc_chisq.h5
MODELFILE_nmo_10kpc_inputdeck.json
```

where `MODELFILE` is the stem of the model file used.

---

## Saved Files

### Ratio Output

The ratio file is saved as:

```text
*_ratio.h5
```

This file contains the coincidence-ratio analysis output produced by:

```python
sf.save_coincidence_output(ratio_result, ratio_outfile)
```

This is the main output used for plotting or analyzing the coincidence-ratio behavior.

---

### Chi-Square and Confidence Contour Output

The chi-square file is saved as:

```text
*_chisq.h5
```

This file contains:

```text
combined_chi_square_grid
ci68
ci99
```

These datasets store the combined chi-square grid and the corresponding 68% and 99% confidence contours.

The file is written by:

```python
save_combined_chi_square(
    chisq_outfile,
    combinedChiSquare,
    ci68Combined,
    ci99Combined
)
```

This is the output used for plotting the confidence contours at each distance.

---

### Input Deck

The input deck is saved as:

```text
*_inputdeck.json
```

This file records the configuration used for the run, including:

```text
modelName
modelPath
modelFileUsed
minMultiplicity
maxMultiplicity
massHierarchy
nGen
distance
distanceLabel
nModules
hitFilePath
spectraFilePath
bkgFilePath
baseFolderPath
```

This makes each run reproducible.

---

## Running Multiple Distances

The `--distance` argument accepts more than one value.

For example:

```bash
--distance 3.0 5.0 10.0 15.0 20.0 25.0
```

This runs the full calculation separately for:

```text
3 kpc
5 kpc
10 kpc
15 kpc
20 kpc
25 kpc
```

Each distance gets its own output folder.

---

## Running Multiple Oscillation Scenarios

To run all oscillation scenarios, use:

```bash
--massHierarchy all
```

This runs:

```text
nmo
imo
no_osc
```

The output will be separated into different folders:

```text
output/
└── modelName/
    ├── nmo/
    ├── imo/
    └── no_osc/
```

---

## Running One Model File vs. a Model Directory

The `--modelPath` argument can point to either a single model file or a directory.

### Single file

```bash
--modelPath **/path/to/model_file.h5**
```

The script runs only that file.

### Directory

```bash
--modelPath **/path/to/model_directory/**
```

The script finds all `.h5` files in the directory and runs them one by one.

---

## Minimal Command Template

Use this template if running manually without the Bash script:

```bash
python3 snecc_main_optimized.py \
  --modelName **MODEL_NAME** \
  --modelPath **/path/to/model_file_or_directory/** \
  --minMultiplicity **2** \
  --maxMultiplicity **15** \
  --massHierarchy **nmo** \
  --nGen **5255019** \
  --distance **3.0 5.0 10.0** \
  --nModules **15000** \
  --hitFilePath **/path/to/hit_file.csv** \
  --spectraFilePath **/path/to/spectra_grid.h5** \
  --output_dir **output/**
```

Replace every bold value with the correct value for your system.

---

## Summary

This pipeline runs the IceCube supernova coincidence-ratio energy-constraint analysis for one or more SNEWPY model files. It supports multiple distances and mass-hierarchy options, saves ratio plot data, saves combined chi-square confidence contour data, and writes an input deck for each run so the analysis can be reproduced later.
