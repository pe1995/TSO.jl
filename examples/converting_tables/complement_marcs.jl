using Pkg; Pkg.activate("."); 
using TSO
using PyPlot
using Glob
using Serialization

## Expected duration for normal tables: ~20 min for square, ~20 min for complement
begin
    ## Load an equation of state that provides the internal energy
    eos = @axed reload(SqEoS, abspath("../../../tests/TSO_sun_Magg_v10.2/eos.hdf5"))

    paths = glob("OS_table*", "OPAC-for-3D/april23/MaggZ0.0a0.0/")  ## Paths to the Tables
    #paths = glob("OS_table*", "../../../create_tables/MARCS/OPAC-for-3D/Z0.0a0.0")  ## Paths to the Tables
    paths = glob("OS_table*", "../../../create_tables/MARCS/OPAC-for-3D/MaggZ0.0a0.0")  ## Paths to the Tables

    mos   = MARCSOpacity(paths...)                     ## Read the raw tables
    m_int = uniform(mos...)                            ## Interpolate to square T-rho grid

    ## Interpolate the tables to new eos
    neweos, newopa, newopa_c, newopa_l, newopa_s = complement(m_int, eos, unify=false)

    TSO.set_limits!(@axed(neweos), newopa)
    TSO.set_limits!(@axed(neweos), newopa_c)
    TSO.set_limits!(@axed(neweos), newopa_l)


    ## Save everything in the usual TSO format
    dname = "TSO_MARCS_v1.1"

    save(neweos,   "combined_eos_marcs.hdf5")
    save(newopa,   "combined_opacities_marcs.hdf5")
    save(newopa_c, "combined_Copacities_marcs.hdf5")
    save(newopa_l, "combined_Lopacities_marcs.hdf5")
    save(newopa_s, "combined_Sopacities_marcs.hdf5")

    ## move
    !isdir(dname) && mkdir(dname)
    mv("combined_eos_marcs.hdf5",        joinpath(dname, "combined_eos.hdf5"),        force=true)
    mv("combined_opacities_marcs.hdf5",  joinpath(dname, "combined_opacities.hdf5"),  force=true)
    mv("combined_Copacities_marcs.hdf5", joinpath(dname, "combined_Copacities.hdf5"), force=true)
    mv("combined_Lopacities_marcs.hdf5", joinpath(dname, "combined_Lopacities.hdf5"), force=true)
    mv("combined_Sopacities_marcs.hdf5", joinpath(dname, "combined_Sopacities.hdf5"), force=true)
end


#= 
v0.5: Use the EoS that was available on gemini at the moment. (TSO_sun_Magg_v10.2)
v1.X: New opacities for Magg, with old EoS (TSO_sun_Magg_v10.2)
=#