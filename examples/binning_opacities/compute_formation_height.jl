using Pkg; Pkg.activate(".")
using TSO
using DelimitedFiles


TSO.load_TS()
TSO.load_wrapper()



# The EoS has already been smoothed in the running process
table_folder = joinpath("tables/TSO_AESOPUS_v1.1")
opacities    = reload(SqOpacity, joinpath(table_folder, "combined_opacities.hdf5"))
eos          = reload(SqEoS,     joinpath(table_folder, "combined_eos.hdf5"))



# For any model, compute the internal energy given T, rho by bisecting using the EoS
solar_model        = reverse(readdlm("stagger_av.dat", skipstart=0), dims=1)
solar_model[:, 1] .= -solar_model[:, 1]
solar_model[:, 3] .= exp.(solar_model[:, 3])
z, lnρ, lnT        = solar_model[:, 1], log.(solar_model[:, 3]), log.(solar_model[:, 2]) 
model              = TSO.lnT_to_lnEi(solar_model, eos)


# Recompute the rosseland optical depth scale based on the combined opacities
# This should be more accurate, because the rosseland optical depth requires 
# integration over frequency, which is not available in the individual runs
    #TSO.rosseland_opacity!(eos.lnRoss, eos, opacities; weights=TSO.ω_midpoint(opacities.λ))
    #TSO.save(eos, joinpath(table_folder, "combined_ross_eos.hdf5"))


# Compute the optical depth scale + the formation height
τ_ross, τ_λ = optical_depth(eos, opacities, model)
d_ross, d_κ = formation_height(model, eos, opacities, τ_ross, τ_λ)



# Save the results
formation_opacities = SqOpacity(d_κ, d_ross, opacities.src, opacities.λ, true)
save(formation_opacities, joinpath(table_folder, "combined_formation_opacities.hdf5"))