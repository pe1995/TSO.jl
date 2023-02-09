struct Model1D{T} <:AbstractModel
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
    logg ::T
end



#= Constructors =#

Model1D(eos, z, lnρ, lnT; logg=log10(2.75e4)) = Model1D(z, lnρ, lnT, lnT_to_lnEi(eos, lnρ, lnT), logg)
Average3D(eos, path; logg=log10(2.75e4))      = begin
    model = Average3D(path)
    e     = lnT_to_lnEi(eos, model.lnρ, model.lnT)
    Model1D(model.z, model.lnρ, model.lnT, e, logg)
end
Average3D(path; logg=log10(2.75e4)) = begin
    solar_model        = reverse(readdlm(path, skipstart=0), dims=1)
    solar_model[:, 1] .= -solar_model[:, 1]
    solar_model[:, 3] .= exp.(solar_model[:, 3])
    z, lnρ, lnT        = solar_model[:, 1], log.(solar_model[:, 3]), log.(solar_model[:, 2]) 
    Model1D(z, lnρ, lnT, zeros(length(z)), logg)
end




#= Model conversions =#

"""
Convert z,T,ρ to z,E,ρ model based on EoS.
"""
lnT_to_lnEi(eos::AxedEoS, lnρ, lnT) = lookup(eos, :lnEi, lnρ, lnT)
lnT_to_lnEi(eos::RegularEoSTable, lnρ, lnT) = lookup(AxedEoS(eos), :lnEi, lnρ, lnT)
