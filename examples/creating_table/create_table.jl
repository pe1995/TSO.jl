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


TSO.load_TS()           # Set the TS path
TSO.load_wrapper()      # Set the Wrapper path
TSO.move_output()       # move previous output to garbage folder in wrapper
TSO.import_wrapper()    # import the python modules of the wrapper, make them avail. to TSO.jl


time = if length(ARGS) >= 1
    ARGS[1]
else
    nothing
end

@show time

##################################################################################################

# Creating initial set of tables
lnT = range(log(1.1e3), log(5e5); length=191)
lnρ = range(log(1e-15), log(1e-3); length=191)
TSO.write_as_stagger(Float64[lnT...], Float64[lnρ...])


# Import the setup class from the wrapper
source_folder = TSO.@inWrapper "source"
TSO.@pythonHelp compute_opac source_folder 


# Basic setup
begin
    ## config file used to setup TurboSpectrum (TS) run
    debug = 1

    ## TS root directory
    ts_root = TSO.@inTS ""
   
    ## list of model atmospheres to use
    atmos_list = TSO.@inWrapper "example/models/TSO_list.in"

    ## path to model atmospheres 
    atmos_path = TSO.@inWrapper "example/models/"

    ## 'm1d' or 'marcs' -> always Stagger for opacity tables
    atmos_format = "stagger"

    ## path to the linelist(s)
    #linelist = ['/u/peitner/Turbospectrum/TSwrapperOpac/nlte_ges_linelist.txt', '/u/peitner/Turbospectrum/TSwrapperOpac/Hlinedata', '/u/peitner/Turbospectrum/TSwrapperOpac/*GESv5.bsyn']
    #linelist = [ TSO.@inTS("../TSwrapperOpac/nlte_ges_linelist.txt"), TSO.@inTS("../TSwrapperOpac/Hlinedata")]
    linelist = abspath.(["LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list", 
                         "LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list"])
                #"/u/peitner/Turbospectrum/TSwrapperOpac/*GESv5.bsyn"]

    ## Dont use this for opacity tables
    #inputParams_file = 'input_param.txt'
    #nlte_grids_path = './nlteGrids/'
    #nlte_config = [ ]
    #nlte_config = [ H : [ 'H/NLTEgrid_H_MARCS_May-10-2021.bin',  'H/auxData_H_MARCS_May-10-2021.txt, 'H/atom.h20'] ]
    #nlte_config += [ Fe : [ 'Fe/NLTEgrid4TS_Fe_MARCS_May-07-2021.bin','Fe/auxData_Fe_MARCS_May-07-2021.dat', 'Fe/atom.fe607a'] ]
    #nlte_config += [ Ba : [ 'Ba/NLTEgrid_Ba_MARCS_May-10-2021.bin', 'Ba/auxData_Ba_MARCS_May-10-2021.txt', 'Ba/atom.ba111' ] ]

    ## starting wavelenght, AA
    lam_start = 1000
    #lam_start = 26000

    ## last wavelenght, AA
    lam_end =  26000
    #lam_end = 200000

    ## resolution per wavelenght (R capital)
    resolution = 25000
    #resolution = 2500

    setup_input = Dict("debug"         =>debug, 
                        "ts_root"      =>ts_root, 
                        "atmos_path"   =>atmos_path, 
                        "atmos_format" =>atmos_format, 
                        "atmos_list"   =>atmos_list, 
                        "linelist"     =>linelist,  
                        "lam_start"    =>lam_start, 
                        "lam_end"      =>lam_end, 
                        "resolution"   =>resolution)

    ## Create the setup object
    setup       = compute_opac.setup(file=setup_input, mode="MAprovided")
    setup.jobID = "TSO"
    
    wvl_set = "UV_dense"
end

# Execute the code in parallel over all nodes/tasks #################################################
# specified in the slurm environment
# NEW: Directly run the TS as job steps of the current slurm installation
TSO.babsma!(setup, [], timeout=time) 

# Now the table can be read, inverted and run again with opacities
lot = glob("_TSOeos_*_TSO.eos")                                                          # List of EoS tables
eos = TSO.load(TSO.EoSTable, lot)                                                        # Read from file and combine


# Save an intermediate step EoS
TSO.save(eos, "eos_$(wvl_set)_raw_step1.hdf5")


# Apply smoothing if needed
TSO.smooth!(eos)                                                                         # Interpolate where NaN


# Invert the EoS and get new grid
new_eos = TSO.energy_grid(eos)                                                           # Switch from T to E grid
                                                                                         # Could also use energy_grid(), which is simpler

# Save the eos
TSO.save(new_eos, "eos_$(wvl_set)_step1.hdf5")

##################################################################################################
# create the new TS input tables with this eos
ee, rr = TSO.meshgrid(new_eos.lnEi, new_eos.lnRho)
TSO.write_as_stagger(new_eos.lnT, rr)   


# Clean the old variables
eos     = nothing
new_eos = nothing
lot     = nothing


# NEW: Directly run the TS as job steps of the current slurm installation
TSO.babsma!(setup, [], timeout=time) 
TSO.bsyn!(  setup, [], timeout=time)


# Now the table can be read, inverted to save
lot = glob("_TSOeos_*_TSO.eos")                                                          # List of EoS tables
eos, opacities = TSO.load(TSO.EoSTable, TSO.OpacityTable, lot, get_individual=true)      # Return cont/line/scat tables separate


TSO.save(eos, "eos_$(wvl_set)_raw_step2.hdf5")
TSO.save(opacities[1],  "opacities_$(wvl_set)_raw_step2.hdf5")
TSO.save(opacities[2], "Copacities_$(wvl_set)_raw_step2.hdf5")
TSO.save(opacities[3], "Lopacities_$(wvl_set)_raw_step2.hdf5")
TSO.save(opacities[4], "Sopacities_$(wvl_set)_raw_step2.hdf5")


# Apply smoothing if needed
TSO.smooth!(eos, opacities)
TSO.save(eos,          "eos_$(wvl_set)_step2.hdf5")
TSO.save(opacities[1], "opacities_$(wvl_set)_step2.hdf5")


# Cut off the edges and make one uniform Ei grid 
new_eos, new_opacities = TSO.unify(eos, opacities)


# Save the end result
TSO.save(new_eos,          "unified_eos_$(wvl_set)_step2.hdf5")
TSO.save(new_opacities[1], "unified_opacities_$(wvl_set)_step2.hdf5")
TSO.save(new_opacities[2], "unified_Copacities_$(wvl_set)_step2.hdf5")
TSO.save(new_opacities[3], "unified_Lopacities_$(wvl_set)_step2.hdf5")
TSO.save(new_opacities[4], "unified_Sopacities_$(wvl_set)_step2.hdf5")

##################################################################################################
# You can filter the models with this if you did not clean before
#lot = lot[TSO.get_TSO_index.(lot) .<= 100]    # Not needed if you cleaned before
