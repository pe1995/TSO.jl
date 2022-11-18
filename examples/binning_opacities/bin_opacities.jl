using Pkg; Pkg.activate(".")
using TSO


TSO.load_TS()
TSO.load_wrapper()


table_folder = joinpath("tables/sun_Magg")


# TS quantities after post-processing
eos           = TSO.reload(TSO.RegularEoSTable,     joinpath(table_folder, "combined_eos.hdf5"))
opacities     = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_opacities.hdf5"), mmap=false)
opacitiesS    = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_Sopacities.hdf5"), mmap=false)
formOpacities = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_formation_opacities.hdf5"), mmap=false)


# λ Integration weights
weights = TSO.ω_midpoint(opacities.λ)


# Choose the method of Binning. For each bin type, the binning function creates the bins (e.g. bin edges)
# The binning() function then sorts the wavelength points from the opacity table into the correct bins, so that 
# the integrated source function and average opacities can be computed in the box_integrated() function accordingly.
bins_tabgen = TSO.TabgenBinning(TSO.EqualTabgenBins, 
                                    opacities=opacities, formation_opacity=log10.(formOpacities.κ))        # A Tabgen styled binning
bins_stagger = TSO.StaggerBinning(TSO.StaggerBins,                                                         #
                                    opacities=opacities, formation_opacity=-log10.(formOpacities.κ_ross),  #
                                    Nbins=12, λ_low=3.7)                                                   # A Stagger styled, L-shaped binning
bins_co = TSO.Co5boldBinning(TSO.Co5boldBins)                                                              # Fixed Co5bold binning


  
# Sort the wavelength points into the bins based on the chosen bin type
binning = TSO.binning(bins_stagger, opacities, -log10.(formOpacities.κ_ross)) 
#binning = TSO.binning(bins_co, opacities, log10.(formOpacities.κ_ross)) 



# Compute binned quantities
binned_opacities = TSO.box_integrated(binning, weights, eos, opacities, opacitiesS)


# Save everything in the dispatch format
TSO.for_dispatch(eos, binned_opacities)
TSO.save(binned_opacities, "binned_opacities.hdf5")


# Move files to the final folder for dispatch
eos_table_name = "TSO_sun_Magg_v2"
!isdir(eos_table_name) && mkdir(eos_table_name) 
mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)

# Copy the eos for convenience. Usually not a big deal because rather small
cp(joinpath(table_folder, "combined_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true)
