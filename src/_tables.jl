"""
Table containing binned Opacities + Source functions on a regular grid.
"""
struct RegularOpacityTable{T<:AbstractFloat, 
                            NOp, NRoss, NSrc} <: AbstractRegularTable
    κ             ::Array{T, NOp}
    κ_ross        ::Array{T, NRoss}
    src           ::Array{T, NSrc}
    λ             ::Array{T, 1}
    optical_depth ::Bool
end

"""
Table containing EoS properties on a regular grid.
"""
struct RegularEoSTable{F<:AbstractFloat, 
                        NRho, NT, NEi, 
                        NPg, NRoss, NNe} <: AbstractRegularTable
    lnRho  :: Array{F, NRho}
    lnT    :: Array{F, NT}
    lnEi   :: Array{F, NEi}
    lnPg   :: Array{F, NPg}
    lnRoss :: Array{F, NRoss}
    lnNe   :: Array{F, NNe}
end

"""
Table containing binned Opacities + Source functions on an irregular grid.
"""
struct IrregularOpacityTable <: AbstractIrregularTable
end

"""
Table containing EoS properties on an irregular grid.
"""
struct IrregularEoSTable <: AbstractIrregularTable
end

"""
Axis info type.
"""
struct EoSAxis{F<:AbstractFloat, Naxis, I<:Integer}
    values    ::Array{F, Naxis}
    dimension ::I
    length    ::I
    name      ::Symbol
end


#= Aliases =#

OpacityTable = Union{RegularOpacityTable, IrregularOpacityTable} 
EoSTable     = Union{RegularEoSTable, IrregularEoSTable} 
SqOpacity    = RegularOpacityTable
SqEoS        = RegularEoSTable




#= Convenience Constructor functions =#
Opacity(args...; regular=true, kwargs...) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...;     regular=true, kwargs...) = regular ? RegularEoSTable(args...;     kwargs...) : IrregularEoSTable(args...;     kwargs...)

RegularOpacityTable(κ::Array, κ_ross::Array, src::Array; optical_depth=false) = RegularOpacityTable(κ, κ_ross, src, optical_depth) 




#= General functions =#

limits(eos::EoSTable) = begin
    eaxis = ndims(eos.lnT)>1

    var_min = eaxis ? minimum(eos.lnEi) : minimum(eos.lnT)
    var_max = eaxis ? maximum(eos.lnEi) : maximum(eos.lnT)
    
    var_min, var_max, minimum(eos.lnRho), maximum(eos.lnRho)
end

"""
Return Energy axis of the given EoS.
"""
function EnergyAxis(eos::E) where {E<:RegularEoSTable}
    eaxis     = ndims(eos.lnEi) == 1 ? :lnEi : :lnT
    axis_val  = getfield(eos, eaxis)
    axis_dims = ndims(axis_val)
    axis_len  = first(size(axis_val))

    EoSAxis(axis_val, axis_dims, axis_len, eaxis)
end

"""
Return Density axis of the given EoS.
"""
function DensityAxis(eos::E) where {E<:RegularEoSTable}
    axis_val  = getfield(eos, :lnRho)
    axis_dims = ndims(axis_val)
    axis_len  = first(size(axis_val))

    EoSAxis(axis_val, axis_dims, axis_len, :lnRho)
end




#= Lookup functions =#

"""
Lookup "what" in the EoS, return the value spliced as args... indicate.
"""
lookup(eos::EoSTable, what::Symbol, rho::AbstractFloat,  var::AbstractFloat, args...)                           = lookup_function(eos, what, args...)(var, rho)
lookup(eos::EoSTable, what::Symbol, rho::AbstractVector, var::AbstractVector, args...)                          = lookup_function(eos, what, args...).(var, rho)           
lookup(eos::EoSTable, opacities::OpacityTable, what::Symbol, rho::AbstractFloat,  var::AbstractFloat, args...)  = lookup_function(eos, opacities, what, args...)(var, rho)
lookup(eos::EoSTable, opacities::OpacityTable, what::Symbol, rho::AbstractVector, var::AbstractVector, args...) = lookup_function(eos, opacities, what, args...).(var, rho)

lookup_function(eos::E, what::Symbol, args...) where {E<:EoSTable} = begin
    eaxis = ndims(eos.lnEi)==1 
    second_axis = eaxis ? eos.lnEi : eos.lnT

    extrapolate(interpolate((second_axis, eos.lnRho), view(getfield(eos, what), :, :, args...), Gridded(Linear())), Line())
end

lookup_function(eos::E, opacities::O, what::Symbol, args...) where {E<:EoSTable, O<:OpacityTable} = begin
    eaxis = ndims(eos.lnEi)==1 
    second_axis = eaxis ? eos.lnEi : eos.lnT

    ip = extrapolate(interpolate((second_axis, eos.lnRho), log.(view(getfield(opacities, what), :, :, args...)), Gridded(Linear())), Line())
    (x1, x2) -> exp(ip(x1, x2)) 
end

"""
Bisect the EOS to get a value of either of the variables 
on the basis of the other one an one parameter. Enter the names and values as kwargs.
The one that shall be found needs to be given a (start, end).
"""
function bisect(eos::EoSTable, iterations=50, antilinear=false; kwargs...)
    eaxis = ndims(eos.lnT)>1
    bisect_var = eaxis ? :lnEi : :lnT
    given_var  = eaxis ? :lnT  : :lnEi

    # now loop and find the requested values
    bisect_res = 0
    b0, b1     = kwargs[bisect_var]
   
    args = Float64[0.0,0.0]
    args[1] = kwargs[given_var]

    for i in 1:iterations
        bisect_res = (b0+b1)/2
        args[2] = bisect_res

        para_res = lookup(eos, given_var, [kwargs[:lnRho]], [args[2]])
        para_res = size(para_res) == 0 ? para_res : para_res[1]
     
        if !(antilinear)
            if para_res > args[1]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        else
            if para_res < args[1]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        end
    end

    return bisect_res
end

function bisect(lnρ::AbstractVector, lnT::AbstractVector, eos::EoSTable)
    lnEi = similar(lnT)
    emin,emax,_,_ = TSO.limits(eos)

    for i in eachindex(lnρ)
        lnEi[i] = TSO.bisect(eos, lnRho=lnρ[i], lnT=lnT[i], lnEi=[emin, emax])
    end

    #hcat(z, exp.(lnEi), exp.(lnρ));
    lnEi
end





#= Optical depth related functions =#

"""
Compute monochromatic and rosseland optical depth scales of the model 
based on the opacity table
"""
function optical_depth(eos::EoSTable, opacities::OpacityTable, model)
    # Model z, ρ and T in cgs
    z, lnρ, lnT = model[:, 1], log.(model[:, 3]), log.(model[:, 2]) 

    T = eltype(opacities.κ)

    τ_λ    = zeros(T, length(lnρ), length(opacities.λ)) 
    τ_ross = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # For each wavelength we integrate along z, z[1]-> surface, z[end]-> bottom
    for i in eachindex(opacities.λ)
        # Look up the opacity in the table
        κ  .= lookup(eos, opacities, :κ, lnρ, lnT, i)
        ρκ .= exp.(lnρ) .* κ

        # Integrate: τ(z) = [ ∫ ρκ dz ]_z0 ^z
        for j in eachindex(z)
            if j==1 
                τ_λ[1, i] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
            else
                τ_λ[j, i] = τ_λ[j-1, i] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
            end
        end
    end

    # Rosseland optical depth
    κ  .= lookup(eos, opacities, :κ_ross, lnρ, lnT)
    ρκ .= exp.(lnρ) .* κ
    for j in eachindex(z)
        if j==1 
            τ_ross[1] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
        else
            τ_ross[j] = τ_ross[j-1] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
        end
    end

    return τ_ross, τ_λ
end

"""
Compute the formation height + opacity, i.e. the rosseland optical depth
where the monochromatic optical depth is 1.
"""
function formation_height(model, eos::EoSTable, opacities::OpacityTable, τ_ross, τ_λ)
    z, lnρ, lnT = model[:, 1], log.(model[:, 3]), log.(model[:, 2]) 
    
    T = eltype(opacities.κ)

    z   = Base.convert.(T, z)
    lnρ = Base.convert.(T, lnρ)
    lnT = Base.convert.(T, lnT)

    rosseland_depth = zeros(T, size(τ_λ, 2))
    opacity_depth   = zeros(T, size(τ_λ, 2))

    lRoss = log.(τ_ross)
    lλ    = log.(τ_λ)

    t_mono = zeros(T, size(τ_λ, 1))
    r_ross = linear_interpolation(lRoss, lnρ, extrapolation_bc=Line())
    T_ross = linear_interpolation(lRoss, lnT, extrapolation_bc=Line())

    for i in axes(τ_λ, 2)
        t_mono .= lλ[:, i]
        rosseland_depth[i] = exp.(linear_interpolation(t_mono, lRoss, extrapolation_bc=Line())(0.0))
        opacity_depth[i]   = lookup(eos, opacities, :κ, 
                                    [r_ross(rosseland_depth[i])], [T_ross(rosseland_depth[i])], i)[1]
    end

    rosseland_depth, opacity_depth
end


"""
Convert a monochromatic EoS into a binned EoS format
"""
toKappaLog!(opacities, eos) = begin
    opacities.src .= log.(opacities.src)

    for i in eachindex(eos.lnRho)
        opacities.κ[:, i, :] .*= exp(eos.lnRho[i])
    end
    opacities.κ .= log.(opacities.κ) 
end





#= Writing in Dispatch format =#

function tabparam(EOSTableFile, RhoEiRadTableFile, nEiBin, nRhoBin, nRadBins, EiMin, EiMax, RhoMin, RhoMax)
    content="&TABPAR\n"*
            "EOSTableFile = '$(EOSTableFile)'\n"*
            "RhoEiRadTableFile = '$(RhoEiRadTableFile)'\n"*
            "NeTgRadTableFile = 'netg_radtab.dat'\n"*
            "nRhoBin = $(nRhoBin)\n"*
            @sprintf("RhoMin = %.7E\n", RhoMin)*
            @sprintf("RhoMax = %.7E\n", RhoMax)*
            "nEiBin = $(nEiBin)\n"*
            @sprintf("EiMin =  %.7E\n", EiMin)*
            @sprintf("EiMax =  %.7E\n", EiMax)*
            "nRadBins = $(nRadBins)\n"*
            "/"

    open("tabparam.in", "w") do f
        write(f, content)
    end

    nothing
end

tabparam(eos::EoSTable, nradbins; eos_file="eostable.dat", opacity_file="rhoei_radtab.dat") = tabparam(eos_file, opacity_file, size(eos.lnPg)..., nradbins, exp.(TSO.limits(eos))...)

"""
Write the EoS + binned opacities in the same format as in Tabgen, so that it can be read by dispatch.
"""
function for_dispatch(eos::EoSTable, χ, S, ϵ)
    ## Check that the grids really are equally spaced
    d = eos.lnEi[2:end] .- eos.lnEi[1:end-1]
    @assert all(d .≈ first(d)) 

    d = eos.lnRho[2:end] .- eos.lnRho[1:end-1]
    @assert all(d .≈ first(d)) 

    f = FortranFile("rhoei_radtab.dat", "w", access="direct", recl=prod(size(S))*4)
    FortranFiles.write(f, rec=1, ϵ)
    FortranFiles.write(f, rec=2, S)
    FortranFiles.write(f, rec=3, χ)
    close(f)

    eos_table = zeros(eltype(eos.lnT), size(eos.lnT)..., 4)
    eos_table[:, :, 1] = eos.lnPg
    eos_table[:, :, 2] = exp.(eos.lnT)
    eos_table[:, :, 3] = eos.lnNe
    eos_table[:, :, 4] = eos.lnRoss

    f = FortranFile("eostable.dat", "w", access="direct", recl=prod(size(eos.lnT))*4*4)
    FortranFiles.write(f, rec=1, eos_table)
    close(f)

    tabparam(eos, size(χ, 3))
end

for_dispatch(eos::EoSTable, opacities::OpacityTable) = for_dispatch(eos, opacities.κ, opacities.src, opacities.κ_ross)
for_dispatch(eos::EoSTable, opacities::OpacityTable, folder::String) = begin
    for_dispatch(eos, opacities)

    save(opacities, "binned_opacities.hdf5")
    save(eos, "binned_eos.hdf5")

    # Move files to the final folder for dispatch
    eos_table_name = folder
    !isdir(eos_table_name) && mkdir(eos_table_name) 

    mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
    mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
    mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)
    mv("binned_eos.hdf5",       joinpath(eos_table_name, "eos.hdf5"),              force=true);
end