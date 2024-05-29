# On the mixing between flavor singlets in lattice gauge theories coupled to matter fields in multiple representations
This repository contains the code used to prepare the plots and results included in [On the mixing between flavor singlets in lattice gauge theories coupled to matter fields in multiple representations [2405.05765]](https://arxiv.org/abs/2405.05765).

## Instructions: Running the analysis
- Install required dependencies (see below)
- Download the hdf5-file `singlets_smeared.hdf5` from the [Zenodo data release](https://zenodo.org/badge/DOI/10.5281/zenodo.11370542) and place it in `input/hdf5data`
- Download the archive `parameters.zip`, decompress it, and place the directory in `input/`
- Download the file `gradient_flow_results.csv`, and place it in the directory `input/gradient_flow_results/`
- Run the analysis using `bash main.sh` within the top level directory
- The figures and tables can then be found in
    - `output/figures/`
    - `output/tables/`

- If you want to start from the raw logs:
    - Download `logfiles_compressed.zip` and `decompress.sh` from Zenodo
    - Decompress the raw log files by executing `decompress.sh` (Note, that just unzipping the archive is not sufficient since the raw log files are zstd-compressed)
    - Place the decompressed directory in `input`
    - Remove the file `input/hdf5data/singlets_smeared.hdf5` if it exists
    - Download the archive `parameters.zip`, decompress it, and place the directory in `input/`
    - Download the file `gradient_flow_results.csv`, and place it in the directory `input/gradient_flow_results/`
    - Run the analysis using `bash main.sh` within the top level directory

In order to respect the dataset size limit on Zenodo, only the relevant channels (γ5, γ0γ5, γi) are written to the hdf5 file. In order to write all channels to the hdf5 file, set the variable 'write_all_channels_to_hdf5' in the file `MixedRepSinglets/main.jl` to 'true'.

## Run times (rough estimate)

- Running from the hdf5 files takes around 20 minutes on a laptop with an i7-1355U CPU
- Running from the raw log files takes around 90 minutes on a laptop with an i7-1355U CPU

## Warning

The code in this repository has only been tested on the specific dataset provided here. It is not intended to be easily generalizable to arbitrary datasets. The analysis parameters are hard-coded in the directory `input/parameters/`.

## Plots

The plots are made using [Plots.jl](https://zenodo.org/record/7994271) via the [PGFPlotsX](https://github.com/KristofferC/PGFPlotsX.jl) backend which requires a LaTeX installation with the PGFPlots package.

## Requirements

- Python 3.10 (see `requirements.txt` for the required packages)
- julia 1.10
- LaTeX (including PGFPlots)
- zstd (for decompressing the raw log files)

### Installing requirements using Conda

All requirements apart from LaTeX may be installed from Conda by running

    conda env create -f environment.yml

On macOS on Apple silicon processors,
it may be necessary to have Rosetta enabled,
and to specify to use x86-64 packages,
by running

    conda env create --platform osx-64 -f environment.yml

Before running the analysis,
the environment should be activated using

    conda activate Sp2N_mixed_representation_singlets_2405.05765

### Reproducibility and reusability

The correlator fitting of this analysis makes use of `scipy.optimize.curve_fit`; as such, results are not bitwise reproducible between different CPU architectures.
