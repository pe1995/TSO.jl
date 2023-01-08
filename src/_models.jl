struct Model1D{T} <:AbstractModel
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
end



#= Constructors =#

Model1D(eos, z, lnρ, lnT) = Model1D(z, lnρ, lnT, lnT_to_lnEi(eos, lnρ, lnT))
Average3D(eos, path) = begin
    solar_model        = reverse(readdlm("stagger_av.dat", skipstart=0), dims=1)
    solar_model[:, 1] .= -solar_model[:, 1]
    solar_model[:, 3] .= exp.(solar_model[:, 3])
    z, lnρ, lnT        = solar_model[:, 1], log.(solar_model[:, 3]), log.(solar_model[:, 2]) 
    
    Model1D(eos, z, lnρ, lnT)
end




#= Model conversions =#

"""
Convert z,T,ρ to z,E,ρ model based on EoS.
"""
lnT_to_lnEi(eos, lnρ, lnT) = lookup(eos, :lnEi, lnρ, lnT)
