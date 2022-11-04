using Pkg; Pkg.activate(".")
using TSO


TSO.load_TS()
TSO.load_wrapper()


# TS quantities after post-processing
eos           = TSO.reload(TSO.RegularEoSTable,     "unified_eos_step2.hdf5")
opacities     = TSO.reload(TSO.RegularOpacityTable, "unified_opacities_step2.hdf5")
opacitiesS    = TSO.reload(TSO.RegularOpacityTable, "unified_Sopacities_step2.hdf5")
formOpacities = TSO.reload(TSO.RegularOpacityTable, "unified_formation_opacities_step2.hdf5")


# λ Integration weights
weights = TSO.ω_midpoint(opacities.λ)


bins_tabgen  = TSO.TabgenBinning(TSO.EqualTabgenBins, opacities=opacities, formation_opacity=log10.(formOpacities.κ))          # A Tabgen styled binning
bins_stagger = TSO.StaggerBinning(TSO.StaggerBins,    opacities=opacities, formation_opacity=log10.(formOpacities.κ), Nbins=8) # A Stagger styled, L-shaped binning
binning      = TSO.binning(bins_stagger, opacities, log10.(formOpacities.κ))                                                   # Sort the wavelength points into the bins based on the chosen bin type


# Compute binned quantities
binned_opacities = TSO.box_integrated(binning, weights, eos, opacities, opacitiesS)


# Save everything in the dispatch format
TSO.for_dispatch(eos, binned_opacities)

eos_table_name = "TSO_table_v1"
!isdir(eos_table_name) && mkdir(eos_table_name) 
mv("tabparam.in",      joinpath(eos_table_name,"tabparam.in"),      force=true)
mv("eostable.dat",     joinpath(eos_table_name,"eostable.dat"),     force=true)
mv("rhoei_radtab.dat", joinpath(eos_table_name,"rhoei_radtab.dat"), force=true);
