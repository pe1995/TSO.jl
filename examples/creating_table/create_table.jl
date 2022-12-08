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
lnT = range(log(1.1e3), log(5.5e5); length=200)
lnρ = range(log(1e-15), log(1e-3); length=200)
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
    linelist = abspath.([#"LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list", 
                         "LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
                         #"LINE-LISTS/25500-200000_cut-4/atom_25500-200000.list",
                         "LINE-LISTS/ADDITIONAL-LISTS/Hlinedata"])
                #"/u/peitner/Turbospectrum/TSwrapperOpac/*GESv5.bsyn"]

    ## Dont use this for opacity tables
    #inputParams_file = 'input_param.txt'
    #nlte_grids_path = './nlteGrids/'
    #nlte_config = [ ]
    #nlte_config = [ H : [ 'H/NLTEgrid_H_MARCS_May-10-2021.bin',  'H/auxData_H_MARCS_May-10-2021.txt, 'H/atom.h20'] ]
    #nlte_config += [ Fe : [ 'Fe/NLTEgrid4TS_Fe_MARCS_May-07-2021.bin','Fe/auxData_Fe_MARCS_May-07-2021.dat', 'Fe/atom.fe607a'] ]
    #nlte_config += [ Ba : [ 'Ba/NLTEgrid_Ba_MARCS_May-10-2021.bin', 'Ba/auxData_Ba_MARCS_May-10-2021.txt', 'Ba/atom.ba111' ] ]

    ## starting wavelenght, AA
    #lam_start = 1000
    lam_start = 2500

    ## last wavelenght, AA
    #lam_end =  26000
    lam_end = 4000

    ## resolution per wavelenght (R capital)
    resolution = 325000
    #resolution = 2500

    tmolim = 20000.0

    @info "Chosen λ step + Number of points: $(TSO.ΔΛ(lam_start,lam_end,resolution)), $(TSO.N_Λ(lam_start,lam_end,resolution))"

    setup_input = Dict("debug"         =>debug, 
                        "ts_root"      =>ts_root, 
                        "atmos_path"   =>atmos_path, 
                        "atmos_format" =>atmos_format, 
                        "atmos_list"   =>atmos_list, 
                        "linelist"     =>linelist,  
                        "lam_start"    =>lam_start, 
                        "lam_end"      =>lam_end, 
                        "resolution"   =>resolution,
                        "TMOLIM"       =>tmolim)

    ## Create the setup object
    setup       = compute_opac.setup(file=setup_input, mode="MAprovided")
    setup.jobID = "TSO"
    
    wvl_set = "UVD_Magg_v1.1"

    magg_2022 = [(TSO.atom_id(:H ), 12.0),
                 (TSO.atom_id(:C ), 8.56),
                 (TSO.atom_id(:N ), 7.98),
                 (TSO.atom_id(:O ), 8.77),
                 (TSO.atom_id(:F ), 4.40),
                 (TSO.atom_id(:Ne), 8.15),
                 (TSO.atom_id(:Na), 6.29),
                 (TSO.atom_id(:Mg), 7.55),
                 (TSO.atom_id(:Al), 6.43),
                 (TSO.atom_id(:Si), 7.59),
                 (TSO.atom_id(:P ), 5.41),
                 (TSO.atom_id(:S ), 7.16),
                 (TSO.atom_id(:Cl), 5.25),
                 (TSO.atom_id(:Ar), 6.50),
                 (TSO.atom_id(:K ), 5.14),
                 (TSO.atom_id(:Ca), 6.37),
                 (TSO.atom_id(:Sc), 3.07),
                 (TSO.atom_id(:Ti), 4.94),
                 (TSO.atom_id(:V ), 3.89),
                 (TSO.atom_id(:Cr), 5.74),
                 (TSO.atom_id(:Mn), 5.52),
                 (TSO.atom_id(:Fe), 7.50),
                 (TSO.atom_id(:Co), 4.95),
                 (TSO.atom_id(:Ni), 6.24)]
    
    # Pick the abundances
    abundances = magg_2022
    @info "Modify the following abundances: (Species, ID, Abundance)"
    for a in abundances
        @info "$(TSO.id_atom(a[1])) - $(a[1]) - $(a[2])"
    end
end
# Execute the code in parallel over all nodes/tasks #################################################
# specified in the slurm environment
# NEW: Directly run the TS as job steps of the current slurm installation
TSO.babsma!(setup, abundances, timeout=time) 

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
TSO.babsma!(setup, abundances, timeout=time) 
TSO.bsyn!(  setup, abundances, timeout=time)


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
