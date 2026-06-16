#!/usr/bin/env python3

import os
import sys
import json
import argparse
import traceback
import h5py
from pathlib import Path

import snecc_functions as sf


# ============================================================
# Utilities
# ============================================================

def fatal_error(stage, err):
    print("\n" + "=" * 70)
    print(f"❌ FATAL ERROR during: {stage}")
    print("=" * 70)
    print(str(err))
    print("\nFull traceback:\n")
    traceback.print_exc()
    print("=" * 70)
    sys.exit(1)


def resolve_mass_hierarchy(option):
    if option == "all":
        return ["nmo", "imo", "no_osc"]
    return [option]


def format_distance_label(distance):
    """
    Format distance for folder/file names.

    This preserves clean labels:
      3.0  -> 3kpc
      3.5  -> 3.5kpc
      10.0 -> 10kpc
    """
    if float(distance).is_integer():
        return f"{int(distance)}kpc"
    return f"{distance:g}kpc"


def save_combined_chi_square(outfile_path, grid, ci68, ci99):
    with h5py.File(outfile_path, "w") as f:
        f.create_dataset("combined_chi_square_grid", data=grid)
        f.create_dataset("ci68", data=ci68)
        f.create_dataset("ci99", data=ci99)


# ============================================================
# Core Runner
# ============================================================

def run_single_model_file(model_file, args, mass_hierarchy, distance):

    distance_label = format_distance_label(distance)

    print("\n" + "-" * 60)
    print(f"Model file: {model_file}")
    print(f"Oscillation: {mass_hierarchy}")
    print(f"Distance: {distance_label}")
    print("-" * 60)

    try:

        # Directory structure uses USER modelName
        distance_folder = os.path.join(
            args.output_dir,
            args.modelName,
            mass_hierarchy,
            distance_label
        )
        os.makedirs(distance_folder, exist_ok=True)

        # Output base name uses file stem
        file_stem = Path(model_file).stem
        base_filename = f"{file_stem}_{mass_hierarchy}_{distance_label}"

        ratio_outfile = os.path.join(distance_folder, base_filename + "_ratio.h5")
        chisq_outfile = os.path.join(distance_folder, base_filename + "_chisq.h5")
        inputdeck_outfile = os.path.join(distance_folder, base_filename + "_inputdeck.json")

        # Load events
        print("Loading events...")
        events = sf.load_events(args.hitFilePath)

        # Build distribution
        print("Building SNEWPY distribution...")
        snewpy_distribution = sf.build_distribution(
            model_file,
            mass_hierarchy=mass_hierarchy
        )

        # Ratio calculation
        print("Calculating coincidence ratios...")
        ratio_result = sf.calculate_coincidence_ratio(
            events=events,
            snewpy_distribution=snewpy_distribution,
            multiplicity_max=args.maxMultiplicity,
            multiplicity_min=args.minMultiplicity,
            filename=args.spectraFilePath,
            Ngen=args.nGen,
            distance=distance,
            N_modules=args.nModules
        )

        print("Saving ratio result...")
        sf.save_coincidence_output(ratio_result, ratio_outfile)

        # Combined chi-square
        print("Computing combined chi-square grid...")
        combinedChiSquare = sf.combined_chi_square_grid(
            ratio_result=ratio_result,
            multiplicity_min=args.minMultiplicity,
            multiplicity_max=args.maxMultiplicity,
            background_file_path=args.bkgFilePath
        )

        print("Computing confidence intervals...")
        chisqCombinedSim, ci68Combined, ci99Combined = (
            sf.calculate_combined_confidence_intervals(
                ratio_result,
                min_multiplicity=args.minMultiplicity,
                max_multiplicity=args.maxMultiplicity
            )
        )

        print("Saving chi-square results...")
        save_combined_chi_square(
            chisq_outfile,
            combinedChiSquare,
            ci68Combined,
            ci99Combined
        )

        # Save input deck
        input_deck = {
            "modelName": args.modelName,
            "modelPath": args.modelPath,
            "modelFileUsed": model_file,
            "minMultiplicity": args.minMultiplicity,
            "maxMultiplicity": args.maxMultiplicity,
            "massHierarchy": mass_hierarchy,
            "nGen": args.nGen,
            "distance": distance,
            "distanceLabel": distance_label,
            "nModules": args.nModules,
            "hitFilePath": args.hitFilePath,
            "spectraFilePath": args.spectraFilePath,
            "bkgFilePath": args.bkgFilePath,
            "baseFolderPath": args.output_dir
        }

        with open(inputdeck_outfile, "w") as f:
            json.dump(input_deck, f, indent=4)

        print("✅ Completed successfully.")

    except Exception as e:
        fatal_error(stage=f"{model_file} | {mass_hierarchy} | {distance_label}", err=e)


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--modelName",
        required=True,
        help="Logical model name, e.g. sukhbold_2015 or tamborra_2014"
    )

    parser.add_argument(
        "--modelPath",
        required=True,
        help="File or directory containing flux models"
    )

    parser.add_argument("--minMultiplicity", type=int, required=True)
    parser.add_argument("--maxMultiplicity", type=int, required=True)

    parser.add_argument(
        "--massHierarchy",
        choices=["nmo", "imo", "no_osc", "all"],
        required=True
    )

    parser.add_argument("--nGen", type=int, required=True)

    # Main change:
    # Accept one or more distances.
    # Examples:
    #   --distance 3
    #   --distance 3 5 10 15 20
    parser.add_argument(
        "--distance",
        type=float,
        nargs="+",
        required=True,
        help="One or more supernova distances in kpc, e.g. --distance 3 5 10"
    )

    parser.add_argument("--nModules", type=int, required=True)
    parser.add_argument("--hitFilePath", required=True)
    parser.add_argument("--spectraFilePath", required=True)
    parser.add_argument("--output_dir", required=True)

    parser.add_argument(
        "--bkgFilePath",
        default="../../OM_bkg_data/30_min_Vessel+PMT_noise/mdom+pmt+30min_bkg_multiplicity.csv"
    )

    args = parser.parse_args()

    try:

        if not os.path.exists(args.modelPath):
            raise FileNotFoundError(f"Model path not found: {args.modelPath}")

        mass_hierarchies = resolve_mass_hierarchy(args.massHierarchy)

        # Directory case
        if os.path.isdir(args.modelPath):

            model_files = [
                os.path.join(args.modelPath, f)
                for f in os.listdir(args.modelPath)
                if f.endswith(".h5")
            ]

            model_files = sorted(model_files)

            if not model_files:
                raise RuntimeError("No .h5 files found in model directory.")

            for model_file in model_files:
                for mh in mass_hierarchies:
                    for distance in args.distance:
                        run_single_model_file(
                            model_file=model_file,
                            args=args,
                            mass_hierarchy=mh,
                            distance=distance
                        )

        # Single file case
        elif os.path.isfile(args.modelPath):

            for mh in mass_hierarchies:
                for distance in args.distance:
                    run_single_model_file(
                        model_file=args.modelPath,
                        args=args,
                        mass_hierarchy=mh,
                        distance=distance
                    )

        else:
            raise ValueError("modelPath must be file or directory.")

    except Exception as e:
        fatal_error(stage="Initialization", err=e)


if __name__ == "__main__":
    main()