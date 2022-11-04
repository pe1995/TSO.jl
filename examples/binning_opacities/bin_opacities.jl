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


# A Tabgen styled binning
bins    = TSO.TabgenBinning(TSO.EqualTabgenBins, opacities=opacities, formation_opacity=log10.(formOpacities.κ))
binning = TSO.binning(bins, opacities, log10.(formOpacities.κ))


# Compute binned quantities
binned_opacities = TSO.box_integrated(binning, weights, eos, opacities, opacitiesS)


# Save everything in the dispatch format
TSO.for_dispatch(eos, binned_opacities)