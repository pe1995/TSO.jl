using Pkg; Pkg.activate(".")
using TSO
using BenchmarkTools


TSO.load_TS()
TSO.load_wrapper()


table_folder = joinpath("tables/sun_Magg_v1.2")


# TS quantities after post-processing
eos           = TSO.reload(TSO.RegularEoSTable,     joinpath(table_folder, "combined_ross_eos.hdf5"))
opacities     = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_opacities.hdf5"), mmap=false)
#opacitiesS    = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_Sopacities.hdf5"), mmap=false)
formOpacities = TSO.reload(TSO.RegularOpacityTable, joinpath(table_folder, "combined_formation_opacities.hdf5"), mmap=true)


# λ Integration weights
weights = TSO.ω_midpoint(opacities.λ)


# Choose the method of Binning. For each bin type, the binning function creates the bins (e.g. bin edges)
# The binning() function then sorts the wavelength points from the opacity table into the correct bins, so that 
# the integrated source function and average opacities can be computed in the box_integrated() function accordingly.
bins_tabgen = TSO.TabgenBinning(TSO.EqualTabgenBins, 
                                    opacities=opacities, formation_opacity=log10.(formOpacities.κ_ross))   # A Tabgen styled binning
bins_stagger = TSO.StaggerBinning(TSO.StaggerBins,                                                         #
                                    opacities=opacities, formation_opacity=-log10.(formOpacities.κ_ross),  #
                                    Nbins=12, λ_low=3.7)                                                   # A Stagger styled, L-shaped binning
bins_co = TSO.Co5boldBinning(TSO.Co5boldBins)                                                              # Fixed Co5bold binning


#= Modifications of the binning =#
## TSO_sun_Magg_v3.1: Use only 8 bins for speed reasons. Also shift the edge a bit lower 
bins_stagger = TSO.StaggerBinning(TSO.StaggerBins,                                                         
                                    opacities=opacities, 
                                    formation_opacity=-log10.(formOpacities.κ_ross),  
                                    Nbins=12, #κ_bins=4,
                                    λ_low=3.6)                 

bins_tabgen = TSO.TabgenBinning(TSO.EqualTabgenBins, 
                                    opacities=opacities, formation_opacity=-log10.(formOpacities.κ_ross), binsize=1.7)   # A Tabgen styled binning
#= End modifications =#
  

# Sort the wavelength points into the bins based on the chosen bin type
binning = TSO.binning(bins_stagger, opacities, -log10.(formOpacities.κ_ross)) 
#binning = TSO.binning(bins_stagger, opacities, -log10.(formOpacities.κ_ross)) 
#binning = TSO.binning(bins_tabgen, opacities, log10.(formOpacities.κ_ross)) 



# Compute binned quantities
binned_opacities = TSO.tabulate(binning, weights, eos, opacities, remove_from_thin=false)

function save(binned_opacities, version)
    # Save everything in the dispatch format
    TSO.for_dispatch(eos, binned_opacities)
    TSO.save(binned_opacities, "binned_opacities.hdf5")


    # Move files to the final folder for dispatch
    eos_table_name = "TSO_sun_Magg_v$(version)"
    !isdir(eos_table_name) && mkdir(eos_table_name) 
    mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
    mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
    mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)

    # Copy the eos for convenience. Usually not a big deal because rather small
    cp(joinpath(table_folder, "combined_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true)
end

function scale!(bo, binned_opacities, what, bins, factor)
    @info "scale by $(log(factor))"
    for bin in bins
        if what == "k"
            bo.κ[:, :, bin] .= log(factor) .+  binned_opacities.κ[:, :, bin]
        elseif what == "s"
            bo.src[:, :, bin] .= log(factor) .+ binned_opacities.src[:, :, bin]
        else
            error("'what' is wrong")
        end
    end
end

save(binned_opacities, "10.2")
#bo = deepcopy(binned_opacities)
#scale!(bo, binned_opacities, "k", [1], 1/2.0)
#scale!(bo, binned_opacities, "k", [4], 2.0)
#save(bo, "9.20")

#=
bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [1], 1/10.0)
scale!(bo, binned_opacities, "s", [4], 10.0)
save(bo, "9.17")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [1], 1/10.0)
scale!(bo, binned_opacities, "k", [4], 10.0)
save(bo, "9.18")
=#
#=
@info binned_opacities.κ[1,1,1]
bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [1,2,3,4], 10.0)
save(bo, "9.2")
@info binned_opacities.κ[1,1,1]

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [1], 10.0)
save(bo, "9.3")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [2], 10.0)
save(bo, "9.4")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [3], 10.0)
save(bo, "9.5")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [4], 10.0)
save(bo, "9.6")



bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [1,2,3,4], 1/10.0)
save(bo, "9.7")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [1], 1/10.0)
save(bo, "9.8")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [2], 1/10.0)
save(bo, "9.9")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [3], 1/10.0)
save(bo, "9.10")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "k", [4], 1/10.0)
save(bo, "9.11")



bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [1,2,3,4], 10.0)
save(bo, "9.12")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [1], 10.0)
save(bo, "9.13")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [2], 10.0)
save(bo, "9.14")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [3], 10.0)
save(bo, "9.15")

bo = deepcopy(binned_opacities)
scale!(bo, binned_opacities, "s", [4], 10.0)
save(bo, "9.16")

=#