using Pkg; Pkg.activate(".")
using TSO
using BenchmarkTools


TSO.load_TS()
TSO.load_wrapper()

name_extension = "DIS_MARCS"
table_folder   = joinpath("tables", "TSO_MARCS_v0.4")


# TS quantities after post-processing
eos           = reload(SqEoS,     joinpath(table_folder, "combined_ross_eos.hdf5")) # ross for others
opacities     = reload(SqOpacity, joinpath(table_folder, "combined_opacities.hdf5"), mmap=true)
formOpacities = reload(SqOpacity, joinpath(table_folder, "combined_formation_opacities.hdf5"), mmap=true)
opacitiesS    = reload(SqOpacity, joinpath(table_folder, "combined_Sopacities.hdf5"), mmap=true)


# λ Integration weights
weights = ω_midpoint(opacities)

# Load a model for the transition between optically thin and thick regime
model = Average3D(eos, "stagger_av.dat")


# Choose the method of Binning. For each bin type, the binning function creates the bins (e.g. bin edges)
# The binning() function then sorts the wavelength points from the opacity table into the correct bins, so that 
# the integrated source function and average opacities can be computed in the box_integrated() function accordingly.
bins_co = Co5boldBinning(TSO.Co5boldBins)                                                                  # Fixed Co5bold binning


#= Modifications of the binning =#
## TSO_sun_Magg_v3.1: Use only 8 bins for speed reasons. Also shift the edge a bit lower 
bins_stagger = StaggerBinning(TSO.StaggerBins,                                                         
                                    opacities=opacities, 
                                    formation_opacity=-log10.(formOpacities.κ_ross),  
                                    Nbins=12, κ_bins=4,
                                    λ_low=3.6)                 

bins_tabgen = TabgenBinning(TSO.EqualTabgenBins, 
                                    opacities=opacities, 
                                    formation_opacity=-log10.(formOpacities.κ_ross), 
                                    binsize=1.7)   # A Tabgen styled binning

bins_beeck = StaggerBinning(TSO.Beeck2012StaggerBins)
bins_muram = MURaMBinning();

#= End modifications =#
  

# Sort the wavelength points into the bins based on the chosen bin type
bin = binning(bins_muram, opacities, -log10.(formOpacities.κ_ross)) 


#for i in 1:12
#    @info "$(count(bin .== i)) wavelengths in bin $(i)"
#end;


# Compute binned quantities
binned_opacities = tabulate(bin, weights, eos, opacities, transition_model=model)



# Save the binned opacities only

function save_table(binned_opacities, version; dispatch=true)
    eos_table_name = "$(name_extension)_v$(version)"
    save(binned_opacities.opacities, "binned_opacities.hdf5")

    !isdir(eos_table_name) && mkdir(eos_table_name) 

    if dispatch
        # Save everything in the dispatch format
        for_dispatch(@axed(eos), binned_opacities)

        # Move files to the final folder for dispatch
        mv("tabparam.in",      joinpath(eos_table_name, "tabparam.in"),      force=true)
        mv("eostable.dat",     joinpath(eos_table_name, "eostable.dat"),     force=true)
        mv("rhoei_radtab.dat", joinpath(eos_table_name, "rhoei_radtab.dat"), force=true)
    end

    # Copy the eos for convenience. Usually not a big deal because rather small
    cp(joinpath(table_folder, "combined_ross_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)
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

save_table(binned_opacities, "0.4.4", dispatch=false)
