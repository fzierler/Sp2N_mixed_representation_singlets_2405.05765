# In order to repsect the dataset size limit on Zenodo, only
# the relevant channels (γ5, γ0γ5, γi) are written to hdf5 files. 
# In order to write all channels to the hdf5 file, set the following
# variable 'write_all_channes_to_hdf5' to 'true'
write_all_channes_to_hdf5 = false

logfiles_path  = normpath("../input/logfiles/")
paramter_path  = normpath("../input/parameters/")
hdf5file_path  = normpath("../input/hdf5data")
output_path    = normpath("../output/")
gradientflowdir= normpath("../input/gradient_flow_results/")
input_h5file   = joinpath(hdf5file_path,"singlets_smeared.hdf5")

start_from_logs    = !isfile(input_h5file)
write_correlator   = true
write_gevp_results = true

if !ispath(paramter_path) || isempty(paramter_path)
    msg = "No parameter files found. Place them in $paramter_path"
    @error msg
    error(msg)
end

if !ispath(hdf5file_path) || !isfile(input_h5file)
    @warn "No hdf5 input file found. I will try and do the analysis starting from the raw logs."

    if !ispath(logfiles_path) || isempty(logfiles_path)
        msg = """ 
        No raw log files, nor hdf5 input file found!
        If you wish to start from the hdf5 file, place 'singlets_smeared.hdf5' in the directory $hdf5file_path
        If you wish to start from the raw log files, place the ceompressed directory within the directory ./input/
        """
        @error msg
        error(msg)
    end

end

include("run_analysis.jl")