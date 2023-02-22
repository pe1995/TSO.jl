#=========================================================== Model Atmosphere structs ==#


struct Model1D{T} <:AbstractModel
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
    logg ::T
end


struct OpticalModel1D{T} <:AbstractModel
    τ    ::Vector{T}
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
    logg ::T
end



#======================================================================= Constructors ==#

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

macro optical(model, eos, opacities)
    m = esc(model)
    e = esc(eos)
    o = esc(opacities)
    quote
        OpticalModel1D(Base.convert(typeof(getfield($m, :z)), rosseland_optical_depth($e, $o, $m)), (getfield($m, p) for p in fieldnames(typeof($m)))...)
    end
end

Model1D(; τ   =nothing,
          z   =nothing,
          lnρ =nothing,
          lnT =nothing,
          lnEi=nothing,
          logg=log10(2.75e4)) = begin
    
    model_arrays = [τ, z, lnρ, lnT, lnEi]
    names        = [:τ, :z, :lnρ, :lnT, :lnEi]

    @assert any(.!isnothing.(model_arrays))

    i_ref = findfirst(.!isnothing.(model_arrays))
    s_ref = size(model_arrays[i_ref])
    t_ref = eltype(model_arrays[i_ref])

    input_arrays = [zeros(t_ref, s_ref...) for _ in names]
    for (i,para) in enumerate(names)
        if isnothing(model_arrays[i])
            continue
        else
            input_arrays[i] = Base.convert.(t_ref, model_arrays[i])
        end
    end

    if isnothing(τ)
        Model1D(input_arrays[2:end]..., logg)
    else
        OpticalModel1D(input_arrays..., logg)
    end
end



#================================================================== Model conversions ==#

"""
Convert z,T,ρ to z,E,ρ model based on EoS.
"""
lnT_to_lnEi(eos::AxedEoS, lnρ, lnT) = lookup(eos, :lnEi, lnρ, lnT)
lnT_to_lnEi(eos::RegularEoSTable, lnρ, lnT) = lookup(AxedEoS(eos), :lnEi, lnρ, lnT)


#=======================================================================================#
