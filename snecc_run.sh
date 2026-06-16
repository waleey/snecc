#!/bin/bash

# ============================================================
# Supernova Energy Constraint Big Run Script
# ============================================================

set -e  # Stop immediately if a command fails
set -o pipefail

echo "============================================================"
echo "Starting SNECC run..."
echo "============================================================"

CMD="python3 snecc_main_optimized.py \
  --modelName sukhbold_sfho \
  --modelPath /Users/walu/icecube/energy_constraints/snecc/data/snewpy_model_data/sukhbold_sfho/ \
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

