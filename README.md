# On the mixing between flavor singlets in lattice gauge theories coupled to matter fields in multiple representations
This repository contains the code used to prepare the plots and results included in [On the mixing between flavor singlets in lattice gauge theories coupled to matter fields in multiple representations [2405.05765]](https://arxiv.org/abs/2405.05765).

## Instructions: Running the analysis
- Install required dependencies (see below)
- Download the hdf5-file `singlets_smeared.hdf5` from the [Zenodo data release]() and place it in `input/hdf5data` (the directory can be modified in `MixedRepSinglets/main.jl`)
- Run the analysis using `bash main.sh` within the top level directory
- The figures and tables can then be found in
    - `output/figures/`
    - `output/tables/`

- If you want to start from the raw logs the variable `start_from_logs` in the file `MixedRepSinglets/main.jl` needs to be set to `true` and the path to the directory containing the decompressed raw logs needs to be provided. Note, that the raw logs are compressed on Zenodo. The variable `start_from_logs` is set to `false` by default.

In order to respect the dataset size limit on Zenodo, only the relevant channels (γ5, γ0γ5, γi) are written to the hdf5 file. In order to write all channels to the hdf5 file, set the variable 'write_all_channels_to_hdf5' in the file `MixedRepSinglets/main.jl` to 'true'.

## Warning

The code in this repository has only been tested on the specific dataset provided here. It is not intended to be easily generalizable to arbitrary datasets. The analysis parameters are hard-coded in the directory `input/parameters/`.

## Plots

The plots are made using [Plots.jl](https://zenodo.org/record/7994271) via the [PGFPlotsX](https://github.com/KristofferC/PGFPlotsX.jl) backend which requires a LaTeX installation with the PGFPlots package.

## Requirements
- Python 3.8 (see `requirements.txt` for the required packages)
- julia 1.10
- LaTeX (including PGFPlots)

