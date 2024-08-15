using Pkg; Pkg.activate("."); Pkg.instantiate();
using MixedRepSinglets
using HiRepParsing
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
using Statistics
using LinearAlgebra
using LsqFit
include("scripts/utils.jl")
include("scripts/write_hdf5.jl")
include("scripts/write_correlatormatrix.jl")
include("scripts/eigenvalues.jl")
include("scripts/massplots.jl")
include("scripts/plotting.jl")
include("scripts/tables.jl")
include("scripts/spectrumplot.jl")
include("scripts/tex_tables.jl")
include("scripts/mixing_angle.jl")
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=1)
Nsmear = collect(0:10:80)

corrfitterpath = joinpath(output_path,"fitresults")
tablepath      = joinpath(output_path,"tables")
tex_tablepath  = joinpath(output_path,"tex_tables")
plotpath       = joinpath(output_path,"plots")

parameterfile      = joinpath(paramter_path,"parameters_smeared.csv")
parameters_fitting = joinpath(paramter_path,"parameters_corrfitter.csv")
parameters_fitting_t0p1 = joinpath(paramter_path,"parameters_corrfitter_t0p1.csv")
parameters_fitting_t0p2 = joinpath(paramter_path,"parameters_corrfitter_t0p2.csv")
parameters_fitting_t0p3 = joinpath(paramter_path,"parameters_corrfitter_t0p3.csv")
parameters_fitting_t0p4 = joinpath(paramter_path,"parameters_corrfitter_t0p4.csv")
parameters_gevp    = joinpath(paramter_path,"parameters_gevp.csv")
gradient_flow_results = joinpath(gradientflowdir,"table1_machine_readable.csv")

ispath(hdf5file_path) || mkpath(hdf5file_path) 

start_from_logs    && main_write_hdf5_logs(Nsmear,logfiles_path,hdf5file_path,parameterfile;filter_channels=!write_all_channes_to_hdf5)
write_correlator   && main_write_correlator_matrices(Nsmear,hdf5file_path)
write_gevp_results && write_eigenvalues(parameters_gevp,hdf5file_path)

function run_corrfitter(parameters_fitting,hdf5file_path,hdf5file,corrfitterpath;resample)
    ispath(corrfitterpath) || mkpath(corrfitterpath)
    resample = resample ? "True" : "False"
    args = `$(abspath(parameters_fitting)) $(abspath(hdf5file_path)) $(abspath(joinpath(hdf5file_path,hdf5file))) $(abspath(corrfitterpath)) $resample`
    try
        run(`python3 scripts/fitting_eigenvalues.py $args`)
    catch
        run(`python  scripts/fitting_eigenvalues.py $args`)
    end
end

run_corrfitter(parameters_fitting,hdf5file_path,"singlets_smeared_eigenvalues.hdf5"     ,corrfitterpath        ;resample=true)
run_corrfitter(parameters_fitting_t0p1,hdf5file_path,"singlets_smeared_eigenvalues_t0p1.hdf5",corrfitterpath*"_t0p1";resample=false)
run_corrfitter(parameters_fitting_t0p2,hdf5file_path,"singlets_smeared_eigenvalues_t0p2.hdf5",corrfitterpath*"_t0p2";resample=false)
run_corrfitter(parameters_fitting_t0p3,hdf5file_path,"singlets_smeared_eigenvalues_t0p3.hdf5",corrfitterpath*"_t0p3";resample=false)
run_corrfitter(parameters_fitting_t0p4,hdf5file_path,"singlets_smeared_eigenvalues_t0p4.hdf5",corrfitterpath*"_t0p4";resample=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,corrfitterpath        ,hdf5file_path,plotpath        ;filename="singlets_smeared_eigenvalues.hdf5",only_singlet=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting_t0p1,corrfitterpath*"_t0p1",hdf5file_path,plotpath*"_t0p1";filename="singlets_smeared_eigenvalues_t0p1.hdf5",only_singlet=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting_t0p2,corrfitterpath*"_t0p2",hdf5file_path,plotpath*"_t0p2";filename="singlets_smeared_eigenvalues_t0p2.hdf5",only_singlet=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting_t0p3,corrfitterpath*"_t0p3",hdf5file_path,plotpath*"_t0p3";filename="singlets_smeared_eigenvalues_t0p1.hdf5",only_singlet=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting_t0p4,corrfitterpath*"_t0p4",hdf5file_path,plotpath*"_t0p4";filename="singlets_smeared_eigenvalues_t0p2.hdf5",only_singlet=false)
write_all_tables(Nsmear,parameters_gevp,parameters_fitting,corrfitterpath,tablepath)
write_all_tables(Nsmear,parameters_gevp,parameters_fitting_t0p1,corrfitterpath*"_t0p1",tablepath*"_t0p1")
write_all_tables(Nsmear,parameters_gevp,parameters_fitting_t0p2,corrfitterpath*"_t0p2",tablepath*"_t0p2")
write_all_tables(Nsmear,parameters_gevp,parameters_fitting_t0p3,corrfitterpath*"_t0p3",tablepath*"_t0p3")
write_all_tables(Nsmear,parameters_gevp,parameters_fitting_t0p4,corrfitterpath*"_t0p4",tablepath*"_t0p4")
write_tex_tables(tablepath,tex_tablepath)
plot_spectrum(tablepath,plotpath,gradient_flow_results)
plot_spectrum(tablepath*"_t0p1",plotpath*"_t0p1",gradient_flow_results)
plot_spectrum(tablepath*"_t0p2",plotpath*"_t0p2",gradient_flow_results)
plot_spectrum(tablepath*"_t0p3",plotpath*"_t0p3",gradient_flow_results)
plot_spectrum(tablepath*"_t0p4",plotpath*"_t0p4",gradient_flow_results)
plot_and_write_mixing_angles(parameters_gevp,hdf5file_path,tablepath,tex_tablepath,plotpath)
write_tex_table_t0_comparison(tablepath,tex_tablepath)