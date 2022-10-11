###########################################################################

"""
Tabgen-styled bins of equal size.
"""
struct EqualTabgenBins{F<:AbstractFloat} <:OpacityBins
    bin_edges::Vector{F}
end

"""
Tabgen-styled bins of almost equal number of points.
"""
struct UniformTabgenBins{F<:AbstractFloat} <:OpacityBins
    bin_edges::Vector{F}
end

struct CustomTabgenBins{F<:AbstractFloat} <:OpacityBins
    bin_edges::Vector{F}
end

EqualTabgenBins()   = EqualTabgenBins(Float64[])
UniformTabgenBins() = UniformTabgenBins(Float64[])

## Alias
TabgenBins = Union{<:EqualTabgenBins, <:UniformTabgenBins, <:CustomTabgenBins}
TabgenBinning(t::Type{<:TabgenBins}, args...; kwargs...) = fill(t, args...; kwargs...)

## Functions for filling the bins
fill(::Type{<:OpacityBins}; kwargs...) = error("Please use specific binning methods.")
fill(::Type{<:CustomTabgenBins}; bin_edges, kwargs...) = CustomTabgenBins(bin_edges)

"""
Equal Tabgen bins based on opacities.
"""
function fill(::Type{<:EqualTabgenBins}; opacities, formation_opacity, Nbins=4, kwargs...)
    # Get the minimum and maximum opacity and decide on a step
    min_opacity  = minimum(formation_opacity)
    max_opacity  = maximum(formation_opacity)

    EqualTabgenBins(collect(range(min_opacity, max_opacity, length=Nbins+1)))
end

"""
Uniformly filled Tabgen bins based on opacities.
"""
function fill(::Type{<:UniformTabgenBins}; opacities, formation_opacity, Nbins=4, kwargs...)
    # sort the opacities and split the array in equal pieces
    sorted_opacities   = sort(formation_opacity)
    splitted_opacities = TSO.split_similar(sorted_opacities, Nbins)
    
    edges = [b[1] for b in splitted_opacities]
    append!(edges, splitted_opacities[end][end])

    UniformTabgenBins(edges)
end

## Functions for computing the binned opacities
"""
Tabgen-styled opacity binning.
    Wavelength points will be sorted into the bins based on their
    opacity at the formation height.
"""
function binning(bins::TabgenBins, opacities, formation_opacity, kwargs...)
    λ = opacities.λ

    # Have the bins been computed already?
    b = if length(bins.bin_edges) == 0
        fill(typeof(bins); opacities=opacities, formation_opacity=formation_opacity, kwargs...)
    else
        bins
    end
    edges = b.bin_edges[1:end-1]

    # the binning (i_bin for every wavelenght)
    binning = zeros(Int, length(λ))
    i_bin ::Union{Nothing,Int}

    for i in eachindex(λ)
        i_bin = findlast(x->(formation_opacity[i]>=x), edges)

        if isnothing(i_bin) # It is smaller than the left edge
            i_bin = 1    
        end

        binning[i] = i_bin
    end

    return binning
end

###########################################################################

"""
Lookup "what" in the EoS, return the value spliced as args... indicate.
"""
function lookup(eos::EoSTable, what::Symbol, rho::AbstractVector, var::AbstractVector, args...)
    eaxis = ndims(eos.lnT)==1 ? false : true
    second_axis = eaxis ? eos.lnEi : eos.lnT

    ip = extrapolate(interpolate((second_axis, eos.lnRho), view(getfield(eos, what), :, :, args...), Gridded(Linear())), Line())
    return ip.(var, rho)
end

"""
Lookup "what" in the EoS, return the value spliced as args... indicate.
"""
function lookup(eos::EoSTable, opacities::OpacityTable, what::Symbol, rho::AbstractVector, var::AbstractVector, args...)
    eaxis = ndims(eos.lnT)==1 ? false : true
    second_axis = eaxis ? eos.lnEi : eos.lnT

    ip = extrapolate(interpolate((second_axis, eos.lnRho), view(getfield(opacities, what), :, :, args...), Gridded(Linear())), Line())
    return ip.(var, rho)
end

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

    hcat(z, exp.(lnEi), exp.(lnρ));
end

### Computation of binned properties ######################################

 

###########################################################################