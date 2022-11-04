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

struct ExactTabgenBins <:OpacityBins
    dbox  ::Float64
    Nbins ::Int
end

struct CustomTabgenBins{F<:AbstractFloat} <:OpacityBins
    bin_edges::Vector{F}
end

struct StaggerBins{F<:AbstractFloat} <:OpacityBins
    bin_edges ::Array{F,2}
end

EqualTabgenBins()   = EqualTabgenBins(Float64[])
UniformTabgenBins() = UniformTabgenBins(Float64[])

## Alias
TabgenBins = Union{<:EqualTabgenBins, <:UniformTabgenBins, <:CustomTabgenBins}
TabgenBinning(t::Type{<:OpacityBins}, args...; kwargs...)  = fill(t, args...; kwargs...)
StaggerBinning(t::Type{<:StaggerBins}, args...; kwargs...) = fill(t, args...; kwargs...)

## Functions for filling the bins
fill(::Type{<:OpacityBins}; kwargs...) = error("Please use specific binning methods.")
fill(::Type{<:CustomTabgenBins}; bin_edges, kwargs...) = CustomTabgenBins(bin_edges)

function fill(::Type{<:StaggerBins}; opacities, formation_opacity, Nbins=6, κ_low=1.5, λ_high=4.0, kwargs...) 
    # L shaped bins of equal size in log λ and log κ
    nbins = Nbins-1
    n_κ_bins = ceil(Int, 3/12 * (nbins))
    n_λ_bins = Int(nbins - n_κ_bins)

    λ_edges = zeros(n_λ_bins+2)
    l = maximum(log10.(opacities.λ)[(formation_opacity .>κ_low ) .& (log10.(opacities.λ) .<λ_high)])
    m = maximum(log10.(opacities.λ))
    r = abs(m - l)
    l -= 0.01*r
    m += 0.01*r
    λ_edges[2:end] .= collect(range(l, m, length=n_λ_bins+1))
    λ_edges[1] = minimum(log10.(opacities.λ))

    # mask of left lambda edge
    mask   = log10.(opacities.λ) .< λ_edges[2]
    k_mask = formation_opacity[mask]
    
    k_edges = zeros(n_κ_bins+1)
    l = minimum(formation_opacity[mask])
    m = maximum(formation_opacity[mask])
    r = abs(m - l)
    l -= 0.01*r
    m += 0.01*r
    k_edges .= collect(range(l, m, length=n_κ_bins+1))

    # Last bin is the one in the top right corner
    bins = zeros(nbins+1 , 4)
    for i in 1:n_λ_bins
        bins[i, 1] = λ_edges[i+1]
        bins[i, 2] = λ_edges[i+2]
        bins[i, 3] = k_edges[1]
        bins[i, 4] = k_edges[2]
    end
    for i in 1:n_κ_bins
        bins[n_λ_bins+i, 1] = λ_edges[1]
        bins[n_λ_bins+i, 2] = λ_edges[2]
        bins[n_λ_bins+i, 3] = k_edges[i]
        bins[n_λ_bins+i, 4] = k_edges[i+1]
    end

    bins[end, 1] = λ_edges[2]
    bins[end, 2] = λ_edges[end]
    bins[end, 3] = k_edges[2]
    bins[end, 4] = k_edges[end]

    StaggerBins(bins)
end

"""
Equal Tabgen bins based on opacities.
"""
function fill(::Type{<:EqualTabgenBins}; opacities, formation_opacity, Nbins=4, binsize=0.0, kwargs...)
    # Get the minimum and maximum opacity and decide on a step
    min_opacity  = minimum(formation_opacity)
    max_opacity  = maximum(formation_opacity)

    return if binsize==0.0
        EqualTabgenBins(collect(range(min_opacity, max_opacity, length=Nbins+1)))
    else
        # Make the last bin bigger
        bins = [min_opacity+binsize*(i-1) for i in 1:Nbins]
        append!(bins, [max_opacity])
        EqualTabgenBins(bins)
    end
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

fill(::Type{<:ExactTabgenBins}; Nbins=4, dbox=0.5, kwargs...) = ExactTabgenBins(dbox, Nbins)

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
    i_bin ::Union{Nothing,Int} = nothing

    for i in eachindex(λ)
        i_bin = findlast(x->(formation_opacity[i]>=x), edges)

        if isnothing(i_bin) # It is smaller than the left edge
            i_bin = 1    
        end

        binning[i] = i_bin
    end

    return binning
end

function binning(bins::ExactTabgenBins, opacities, formation_opacity, kwargs...)
    λ = opacities.λ

    # the binning (i_bin for every wavelenght)
    binning = zeros(Int, length(λ))
    i_bin ::Union{Nothing,Int} = nothing

    for i in eachindex(λ)
        i_bin = floor(Int, 1.5+(formation_opacity[i])/bins.dbox)
        i_bin = min(bins.Nbins, max(1, i_bin))

        binning[i] = i_bin
    end

    return binning
end

function binning(b::StaggerBins, opacities, formation_opacity, kwargs...)
    bins = b.bin_edges
    binning = zeros(Int, size(opacities.λ)...)
    for i in eachindex(opacities.λ)
        for j in axes(bins, 1)
            bin = bins[j, :]
            if (log10(opacities.λ[i]) >= bin[1]) & (log10(opacities.λ[i]) < bin[2])
                if (formation_opacity[i] >= bin[3]) & (formation_opacity[i] < bin[4])
                    binning[i] = j 
                end
            end
        end
    end

    @assert all(binning .!= 0)
    binning
end

"""Compute midpoint-weights from wavelenght array."""
function ω_midpoint(λ) 
    Δλ = zeros(eltype(λ), length(λ)+1)
    Δλ[2:end-1] = (λ[2:end] .- λ[1:end-1]) ./ 2
    Δλ[1]   = Δλ[2]
    Δλ[end] = Δλ[end-1]

    w  = similar(λ)
    for i in eachindex(λ)
        w[i] = Δλ[i] + Δλ[i+1]
    end

    w
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

"""
Loop through eos and integrate quantities in the opacity table according to the chosen binning.
Return ϵ in the κ_ross field of the opacity table.
"""
function box_integrated(binning, weights, eos, opacities, scattering)
    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    EiBins  = length(eos.lnEi)
    λ       = opacities.λ
    T       = eltype(eos.lnRho)

    χBox  = zeros(T, EiBins, rhoBins, radBins)    
    χRBox = zeros(T, EiBins, rhoBins, radBins)
    κBox  = zeros(T, EiBins, rhoBins, radBins)
    SBox  = zeros(T, EiBins, rhoBins, radBins)

    BBox = zeros(T, EiBins, radBins)
    δBBox = zeros(T, EiBins, radBins)

    B  = zeros(T, EiBins, length(λ))
    δB = zeros(T, EiBins, length(λ))

    ρ = exp.(eos.lnRho)

    for i in eachindex(eos.lnRho)

        BBox  .= 0.0
        δBBox .= 0.0

        # compute the Planck function at this density
        B  .= Bν( λ, exp.(eos.lnT[:, i]))
        δB .= δBν(λ, exp.(eos.lnT[:, i]))

        # Integrate the Planck Box
        @inbounds for j in eachindex(λ)
            BBox[:, binning[j]]  .+= weights[j] .*  B[:, j]
            δBBox[:, binning[j]] .+= weights[j] .* δB[:, j]
        end

        # Integrate the other Box quantities
        @inbounds for j in eachindex(λ)
            χBox[:, i, binning[j]]  .+= weights[j] .* opacities.κ[:, i, j]                            .*  B[:, j] ./  BBox[:, binning[j]]
            χRBox[:, i, binning[j]] .+= weights[j] .* 1.0 ./ opacities.κ[:, i, j]                     .* δB[:, j] ./ δBBox[:, binning[j]]
            κBox[:, i, binning[j]]  .+= weights[j] .* (opacities.κ[:, i, j] .- scattering.κ[:, i, j]) .*  B[:, j] ./  BBox[:, binning[j]]
            SBox[:, i, binning[j]]  .+= weights[j] .* B[:, j]
        end

        χBox[:, i, :]  .*= ρ[i]
        κBox[:, i, :]  .*= ρ[i]
        χRBox[:, i, :] ./= ρ[i]
    end

    wthin  = exp.(T(-1.5e7) .* T(1.0) ./χRBox)
    wthick = T(1.0) .- wthin

    # Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    RegularOpacityTable(log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(κBox ./ χBox), log.(SBox), collect(T, 1:radBins), false)
end

###########################################################################
