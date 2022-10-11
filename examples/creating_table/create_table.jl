###########################################################
## Example script for creating an EoS table from regular ##
## grid. Code can be executed on cluster using run.slurm ##
###########################################################

using Pkg; Pkg.activate(".");
using Distributed
using SlurmClusterManager
using LaTeXStrings
using FortranFiles
using StatsBase
using TSO
using DelimitedFiles
using Glob
using Interpolations


TSO.load_TS()       # Set the TS path
TSO.load_wrapper()  # Set the Wrapper path
TSO.move_output()   # move previous output to garbage folder in wrapper


if "SLURM_NTASKS" in keys(ENV)
    addprocs(SlurmManager())
    for i in workers()
        host,pid = fetch(@spawnat i (gethostname(), getpid()))
        @info "Worker $(i) is running on $(host) with ID $(pid)" 
    end
else
    @warn "No Slurm environment detected. Using default addprocs."
    addprocs(3)
end


@everywhere begin
    using Pkg; Pkg.activate("."); 
    using Distributed
    using StatsBase
    using TSO
    using DelimitedFiles
    using Glob
    using Interpolations

    TSO.load_TS()
    TSO.load_wrapper()
end


# Creating initial set of tables
lnT = range(log(1e3), log(5e5); length=100)
lnρ = range(log(1e-15), log(1e-3); length=100)
TSO.write_as_stagger(Float64[lnT...], Float64[lnρ...])


# Import the wrapper everywhere
@everywhere begin 
    source_folder = TSO.@inWrapper "source"
    TSO.@pythonHelp compute_opac source_folder 
end


# Basic setup
begin
    # config file used to setup TurboSpectrum (TS) run
    debug = 1

    # TS root directory
    ts_root = "/u/peitner/Turbospectrum/Turbospectrum_NLTE"
    # list of model atmospheres to use

    # path to model atmospheres (in MARCS .mod format only for now)
    atmos_path = "/u/peitner/Turbospectrum/TurboSpectrum-Wrapper/example/models/"

    # 'm1d' or 'marcs'
    atmos_format = "stagger"
    atmos_list   = "/u/peitner/Turbospectrum/TurboSpectrum-Wrapper/example/models/TSO_list.in"

    # path to the linelist(s)
    #linelist = ['/u/peitner/Turbospectrum/TSwrapperOpac/nlte_ges_linelist.txt', '/u/peitner/Turbospectrum/TSwrapperOpac/Hlinedata', '/u/peitner/Turbospectrum/TSwrapperOpac/*GESv5.bsyn']
    linelist = ["/u/peitner/Turbospectrum/TSwrapperOpac/nlte_ges_linelist.txt", "/u/peitner/Turbospectrum/TSwrapperOpac/Hlinedata"]

    #inputParams_file = 'input_param.txt'
    #nlte_grids_path = './nlteGrids/'
    #nlte_config = [ ]
    #nlte_config = [ H : [ 'H/NLTEgrid_H_MARCS_May-10-2021.bin',  'H/auxData_H_MARCS_May-10-2021.txt, 'H/atom.h20'] ]
    #nlte_config += [ Fe : [ 'Fe/NLTEgrid4TS_Fe_MARCS_May-07-2021.bin','Fe/auxData_Fe_MARCS_May-07-2021.dat', 'Fe/atom.fe607a'] ]
    #nlte_config += [ Ba : [ 'Ba/NLTEgrid_Ba_MARCS_May-10-2021.bin', 'Ba/auxData_Ba_MARCS_May-10-2021.txt', 'Ba/atom.ba111' ] ]

    # starting wavelenght, AA
    lam_start = 3000

    # last wavelenght, AA
    lam_end = 10000

    # resolution per wavelenght (R capital)
    resolution = 1300

    setup_input = Dict("debug"         =>debug, 
                        "ts_root"      =>ts_root, 
                        "atmos_path"   =>atmos_path, 
                        "atmos_format" =>atmos_format, 
                        "atmos_list"   =>atmos_list, 
                        "linelist"     =>linelist,  
                        "lam_start"    =>lam_start, 
                        "lam_end"      =>lam_end, 
                        "resolution"   =>resolution)

    # Create the setup object
    setup       = compute_opac.setup(file=setup_input, mode="MAprovided")
    setup.jobID = "TSO"
end

# Execute the code in parallel over all nodes/tasks
# specified in the slurm environment
# args: setup, atmos number, skip bsyn?
args = [(setup, [i-1], false) for i in 1:length(setup.atmos_list)]
Distributed.pmap(compute_opac.runTSforOpac, args)

#=
# Now the table can be read, inverted and run again with opacities
lot = glob("_TSOeos_*_TSO.eos")               # List of EoS tables
#lot = lot[TSO.get_TSO_index.(lot) .<= 100]    # Not needed if you cleaned before
eos = TSO.read_tables(TSO.EoSTable, list_of_eos_tables)

# Invert the EoS and get new grid
new_eos = TSO.lnT_to_lnEi(d_eos)

# TODO: Apply smoothing if needed
# smooth!(new_eos)

# create the new TS input tables with this eos
ee, rr = TSO.meshgrid(new_eos.lnEi, new_eos.lnRho);
TSO.write_as_stagger(new_eos.lnT, rr)

# Run the code again
args = [(setup, [i-1], false) for i in 1:length(setup.atmos_list)]
Distributed.pmap(compute_opac.runTSforOpac, args)
=#