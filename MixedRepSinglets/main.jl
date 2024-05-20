start_from_logs    = false
write_correlator   = true
write_gevp_results = true

logfiles_path  = normpath("../input/logfiles/")
paramter_path  = normpath("../input/parameters/")
hdf5file_path  = normpath("../input/hdf5data")
output_path    = normpath("../output/")

# In order to repsect the dataset size limit on Zenodo, only_singlet
# the relevant channels (γ5, γ0γ5, γi) are written to hdf5 file sizes. 
# In order to write all channels to the hdf5 file, set the following
# variable 'write_all_channes_to_hdf5' to 'true'
write_all_channes_to_hdf5 = false

include("run_analysis.jl")