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

struct Co5boldBins{F<:AbstractFloat} <:OpacityBins
    bin_edges ::Array{F,2}
end



EqualTabgenBins()   = EqualTabgenBins(Float64[])
UniformTabgenBins() = UniformTabgenBins(Float64[])



## Alias
TabgenBins = Union{<:EqualTabgenBins, <:UniformTabgenBins, <:CustomTabgenBins}
TabgenBinning(t::Type{<:OpacityBins}, args...; kwargs...)  = fill(t, args...; kwargs...)
StaggerBinning(t::Type{<:StaggerBins}, args...; kwargs...) = fill(t, args...; kwargs...)
Co5boldBinning(t::Type{<:Co5boldBins}, args...; kwargs...) = fill(t, args...; kwargs...)
MURaMBinning(args...; kwargs...)                           = fill(CustomTabgenBins, bin_edges=[-99.0,0.0,2.0,4.0,99.0], args...; kwargs...)


## Functions for filling the bins
fill(::Type{<:OpacityBins}; kwargs...) = error("Please use specific binning methods.")
fill(::Type{<:CustomTabgenBins}; bin_edges, kwargs...) = CustomTabgenBins(bin_edges)

fill(::Type{Co5boldBins}; kwargs...) = begin
    bin_edges = zeros(12, 4)

    bin_edges[1,  3] = 0.15;  bin_edges[1,  4] = 99.0;  bin_edges[1,  1] = 0.0;    bin_edges[1,  2] = 5500.0
    bin_edges[2,  3] = 0.15;  bin_edges[2,  4] = 99.0;  bin_edges[2,  1] = 5500.0; bin_edges[2,  2] = 1000000.0
    bin_edges[3,  3] = 0.00;  bin_edges[3,  4] = 0.15;  bin_edges[3,  1] = 0.0;    bin_edges[3,  2] = 6000.0
    bin_edges[4,  3] = 0.00;  bin_edges[4,  4] = 0.15;  bin_edges[4,  1] = 6000.0; bin_edges[4,  2] = 1000000.0
    bin_edges[5,  3] = -0.75; bin_edges[5,  4] = 0.00;  bin_edges[5,  1] = 0.0;    bin_edges[5,  2] = 6500.0
    bin_edges[6,  3] = -0.75; bin_edges[6,  4] = 0.00;  bin_edges[6,  1] = 6500.0; bin_edges[6,  2] = 1000000.0
    bin_edges[7,  3] = -1.50; bin_edges[7,  4] = -0.75; bin_edges[7,  1] = 0.0;    bin_edges[7,  2] = 1000000.0
    bin_edges[8,  3] = -2.25; bin_edges[8,  4] = -1.50; bin_edges[8,  1] = 0.0;    bin_edges[8,  2] = 1000000.0
    bin_edges[9,  3] = -3.00; bin_edges[9,  4] = -2.25; bin_edges[9,  1] = 0.0;    bin_edges[9,  2] = 1000000.0
    bin_edges[10, 3] = -3.75; bin_edges[10, 4] = -3.00; bin_edges[10, 1] = 0.0;    bin_edges[10, 2] = 1000000.0
    bin_edges[11, 3] = -4.50; bin_edges[11, 4] = -3.75; bin_edges[11, 1] = 0.0;    bin_edges[11, 2] = 1000000.0
    bin_edges[12, 3] = -99.0; bin_edges[12, 4] = -4.50; bin_edges[12, 1] = 0.0;    bin_edges[12, 2] = 1000000.0

    Co5boldBins(bin_edges)
end

function fill(::Type{<:StaggerBins}; opacities, formation_opacity, Nbins=6, κ_low=1.5, λ_high=4.0, λ_low=nothing, κ_bins=nothing, kwargs...) 
    # L shaped bins of equal size in log λ and log κ
    nbins = Nbins-1
    n_κ_bins = isnothing(κ_bins) ? ceil(Int, 3/12 * (nbins)) : κ_bins
    n_λ_bins = Int(nbins - n_κ_bins)

    λ_edges = zeros(n_λ_bins+2)
    l = isnothing(λ_low) ? 
            maximum(log10.(opacities.λ)[(formation_opacity .>κ_low ) .& (log10.(opacities.λ) .<λ_high)]) :
            λ_low
    m = maximum(log10.(opacities.λ))
    r = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    λ_edges[2:end] .= collect(range(l, m, length=n_λ_bins+1))
    λ_edges[1] = minimum(log10.(opacities.λ))

    # mask of left lambda edge
    mask   = log10.(opacities.λ) .< λ_edges[2]
    k_mask = formation_opacity[mask]
    
    k_edges = zeros(n_κ_bins+1)
    l = minimum(formation_opacity)
    m = maximum(formation_opacity)
    r = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    l_bin_k_edge = mean(formation_opacity[log10.(opacities.λ) .>= λ_edges[end-1]])

    # Split the part above the first bin equally
    k_splitted  = split_similar(sort(formation_opacity[mask .& (formation_opacity.>l_bin_k_edge)]), n_κ_bins-1)

    k_edges[1]  = l
    k_edges[2]  = l_bin_k_edge
    for b in 2:n_κ_bins
        k_edges[b+1] = last(k_splitted[b-1])
    end
    #k_edges .= collect(range(l, m, length=n_κ_bins+1))
    k_edges[end] = m
    
    # Last bin is the one in the top right corner
    bins = zeros(nbins+1 , 4)
    for i in 1:n_λ_bins
        bins[i, 1] = λ_edges[i+1]
        bins[i, 2] = λ_edges[i+2]
        bins[i, 3] = k_edges[1]
        bins[i, 4] = l_bin_k_edge
    end
    for i in 1:n_κ_bins
        bins[n_λ_bins+i, 1] = λ_edges[1]
        bins[n_λ_bins+i, 2] = λ_edges[2]
        bins[n_λ_bins+i, 3] = k_edges[i]
        bins[n_λ_bins+i, 4] = k_edges[i+1]
    end

    bins[end, 1] = λ_edges[2]
    bins[end, 2] = λ_edges[end]
    bins[end, 3] = l_bin_k_edge
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

function binning(b::B, opacities, formation_opacity, kwargs...) where {B<:Union{<:StaggerBins, <:Co5boldBins}}
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

## Loopup functions
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

lookup(eos::EoSTable, what::Symbol, rho::AbstractFloat, var::AbstractFloat, args...)                          = first(lookup(eos,            what, [rho], [var], args...))
lookup(eos::EoSTable, opacities::OpacityTable, what::Symbol, rho::AbstractFloat, var::AbstractFloat, args...) = first(lookup(eos, opacities, what, [rho], [var], args...))


## Optical depth related functions
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

"""Convert a monochromatic EoS into a binned EoS format"""
toKappaLog!(opacities, eos) = begin
    opacities.src .= log.(opacities.src)

    for i in eachindex(eos.lnRho)
        opacities.κ[:, i, :] .*= exp(eos.lnRho[i])
    end
    opacities.κ .= log.(opacities.κ) 
end



### Computation of binned properties ######################################

## Storage for computation
struct BinningStorage{F<:AbstractFloat}
    wT  :: Vector{F} 
    kT  :: Vector{F} 
    sT  :: Vector{F} 
    BT  :: Vector{F} 
    dBT :: Vector{F} 
    lT  :: Vector{F}
    sizes :: Vector{Int}
end

BinningStorage{T}(Λ) where {T<:AbstractFloat} = begin
    lλ, radBins = size(Λ)
    lΛ = falses(lλ)

    sizes = zeros(Int, radBins)
    for i in 1:radBins
        lΛ .= Λ[:, i]
        sizes[i] = count(lΛ)
    end

    BinningStorage( zeros(T, lλ), 
                    zeros(T, lλ), 
                    zeros(T, lλ), 
                    zeros(T, lλ), 
                    zeros(T, lλ),
                    zeros(T, lλ),
                    sizes)
end 


"""
Compute the intensity for all the density value in the table.
"""
function extrapolate_intensity(model, J, eos)
    lnRho = log.(model[:, 3])
    ip    = linear_interpolation(lnRho, J, extrapolation_bc=Line())
    
    J_new = similar(eos.lnPg)
    for i in eachindex(eos.lnRho)
        J_new[:, i] .= exp(ip(eos.lnRho))
    end

    J_new
end

"""
Loop through eos and integrate quantities in the opacity table according to the chosen binning.
Return ϵ in the κ_ross field of the opacity table.
"""
function box_integrated(binning, weights, eos, opacities, scattering; remove_from_thin=false)
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

    rtest = argmin(abs.(ρ .- test_r))
    ttest = argmin(abs.(eos.lnT[:, rtest] .- test_t))

    # None of the bins is allowed to be empty as this stage
    @assert all(binning .!= 0)

    for i in eachindex(eos.lnRho)

        BBox  .= 0.0
        δBBox .= 0.0

        # compute the Planck function at this density
        B  .= Bν( λ, exp.(eos.lnT[:, i]))
        δB .= δBν(λ, exp.(eos.lnT[:, i]))

        # Integrate the Planck Box
        #@inbounds @views for j in eachindex(λ)
        #    BBox[:, binning[j]]  .+= weights[j] .*  B[:, j]
        #    δBBox[:, binning[j]] .+= weights[j] .* δB[:, j]
        #end
       
        # Integrate the other Box quantities
        @inbounds @views for j in eachindex(λ)
            # Integrate the Planck Box
            BBox[:,  binning[j]] .+= weights[j] .*  B[:, j]
            δBBox[:, binning[j]] .+= weights[j] .* δB[:, j]

            χBox[:, i,  binning[j]] .+= weights[j] .* opacities.κ[:, i, j]                            .*  B[:, j]
            χRBox[:, i, binning[j]] .+= weights[j] .* 1.0 ./ opacities.κ[:, i, j]                     .* δB[:, j]
            κBox[:, i,  binning[j]] .+= weights[j] .* (opacities.κ[:, i, j] .- scattering.κ[:, i, j]) .*  B[:, j] 
            SBox[:, i,  binning[j]] .+= weights[j] .* B[:, j]
        end

        @inbounds @views for j in axes(BBox, 2)
            χBox[:,  i, j] ./= BBox[:, j]
            χRBox[:, i, j] ./= δBBox[:, j]
            κBox[:,  i, j] ./= BBox[:, j]
        end

        if i == rtest
            @info "TestBin: $(test_bin), rho: $(ρ[i])"
            @info "χ: $(log(χBox[ttest, i, test_bin])), χR: $(log(1. /χRBox[ttest, i, test_bin])), simple mean: $(log(mean(χBox[ttest, i, test_bin])))"
            @info "wthin: $(exp(T(-1.5e7) .* T(1.0) ./χRBox[ttest, i, test_bin] *ρ[i]))"
            #@info "BBox: $(BBox[ttest, test_bin]), two B point: $(weights[findfirst(Λ[:,test_bin])] .*B[ttest, findfirst(Λ[:,test_bin])]), $(weights[findlast(Λ[:,test_bin])] .*B[ttest, findlast(Λ[:,test_bin])])"
        end

        χBox[:, i, :]  .*= ρ[i]
        κBox[:, i, :]  .*= ρ[i]
        χRBox[:, i, :] ./= ρ[i]
    end

    wthin  = exp.(T(-1.5e7) .* T(1.0) ./χRBox)
    #wthin  = exp.(T(-2) .* T(1.0) ./χRBox)
    wthick = T(1.0) .- wthin

    # Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    # Opacity
    opacity_table = remove_from_thin ? 
                    log.(wthin .* κBox .+ wthick .* (T(1.0) ./ χRBox)) : 
                    log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox))

    # ϵ table
    ϵ_table = log.(κBox ./ χBox)

    # Source function table --> S                = x/x+s * B + s/x+s * J 
    #                       --> thermal emission = x/x+s * B
    S_table = log.(SBox)
    #S_table = log.(κBox ./ χBox .* SBox)

    RegularOpacityTable(opacity_table, ϵ_table, S_table, collect(T, 1:radBins), false)
end

"""
Loop through eos and integrate quantities in the opacity table according to the chosen binning.
Return ϵ in the κ_ross field of the opacity table.
"""
function box_integrated_v2(binning, weights, eos, opacities, scattering; remove_from_thin=false)
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

    rtest = argmin(abs.(ρ .- test_r))
    ttest = argmin(abs.(eos.lnT[:, rtest] .- test_t))

    # None of the bins is allowed to be empty as this stage
    @assert all(binning .!= 0)

    Λ = falses(length(λ), radBins)
    for i in 1:radBins
        Λ[:, i] .= binning .== i
    end
    lΛ = falses(length(λ))

    ### Storage for the loop variables to avoid allocations
    store = BinningStorage{T}(Λ)

    for i in eachindex(eos.lnRho)

        BBox  .= 0.0
        δBBox .= 0.0

        # compute the Planck function at this density
        #B  .= Bν( λ, exp.(eos.lnT[:, i]))
        #δB .= δBν(λ, exp.(eos.lnT[:, i]))

        # Integrate the Planck Box
        #=@inbounds for j in 1:radBins
            copyto!(lΛ, view(Λ, :, j))
            @inbounds for k in eachindex(eos.lnEi)
                copyto!(store.wT,  view(weights ,lΛ))
                copyto!(store.BT,  view(B,  k, lΛ))
                copyto!(store.dBT, view(δB, k, lΛ))

                @inbounds for l in 1:store.sizes[j]
                    BBox[ k, j] += store.wT[l] .* store.BT[l]
                    δBBox[k, j] += store.wT[l] .* store.dBT[l]
                end
            end
        end=#
       
        # Integrate the other Box quantities
        @inbounds for j in 1:radBins
            copyto!(lΛ, view(Λ, :, j))
            @inbounds for k in eachindex(eos.lnEi)
                clear!(store)
                copyto!(store.wT,  view(weights, lΛ))
                copyto!(store.lT,  view(opactities.λ, lΛ))
                #copyto!(store.BT,  view(B,  k, lΛ))
                #copyto!(store.dBT, view(δB, k, lΛ))
                copyto!(store.kT,  view(opacities.κ,  k, i, lΛ))
                copyto!(store.sT,  view(scattering.κ, k, i, lΛ))

                @inbounds for l in 1:store.sizes[j]
                    χBox[ k, i, j]  += store.wT[l] .* store.kT[l] .*                  Bν( store.lT[l], T[k, i]) #BT[l]  
                    χRBox[k, i, j]  += store.wT[l] .* 1.0 ./ store.kT[l]           .* δBν(store.lT[l], T[k, i]) #dBT[l] 
                    κBox[ k, i, j]  += store.wT[l] .* (store.kT[l] .- store.sT[l]) .* Bν( store.lT[l], T[k, i]) #BT[l]  
                    SBox[ k, i, j]  += store.wT[l] .*                                 Bν( store.lT[l], T[k, i]) #BT[l] 
                    BBox[ k, j]     += store.wT[l] .*                                 Bν( store.lT[l], T[k, i]) #BT[l]
                    δBBox[k, j]     += store.wT[l] .*                                 δBν(store.lT[l], T[k, i]) #dBT[l]
                end
                χBox[k, i, j]  /=   BBox[k, j]
                κBox[k, i, j]  /=   BBox[k, j]
                χRBox[k, i, j] /=  δBBox[k, j]
            end
        end

        if i == rtest
            @info "TestBin: $(test_bin), rho: $(ρ[i])"
            @info "χ: $(log(χBox[ttest, i, test_bin])), χR: $(log(1. /χRBox[ttest, i, test_bin])), simple mean: $(log(mean(χBox[ttest, i, test_bin])))"
            @info "wthin: $(exp(T(-1.5e7) .* T(1.0) ./χRBox[ttest, i, test_bin] *ρ[i]))"
            @info "BBox: $(BBox[ttest, test_bin]), two B point: $(weights[findfirst(Λ[:,test_bin])] .*B[ttest, findfirst(Λ[:,test_bin])]), $(weights[findlast(Λ[:,test_bin])] .*B[ttest, findlast(Λ[:,test_bin])])"
        end

        χBox[:, i, :]  .*= ρ[i]
        κBox[:, i, :]  .*= ρ[i]
        χRBox[:, i, :] ./= ρ[i]
    end

    wthin  = exp.(T(-1.5e7) .* T(1.0) ./χRBox)
    #wthin  = exp.(T(-2) .* T(1.0) ./χRBox)
    wthick = T(1.0) .- wthin

    # Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    # Opacity
    opacity_table = remove_from_thin ? 
                    log.(wthin .* κBox .+ wthick .* (T(1.0) ./ χRBox)) : 
                    log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox))

    # ϵ table
    ϵ_table = log.(κBox ./ χBox)

    # Source function table --> S                = x/x+s * B + s/x+s * J 
    #                       --> thermal emission = x/x+s * B
    S_table = log.(SBox)
    #S_table = log.(κBox ./ χBox .* SBox)

    RegularOpacityTable(opacity_table, ϵ_table, S_table, collect(T, 1:radBins), false)
end

function box_integrated_v3(binning, weights, eos, opacities, scattering; remove_from_thin=false)
    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    EiBins  = length(eos.lnEi)
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, EiBins, rhoBins, radBins)    
    χRBox  = zeros(T, EiBins, rhoBins, radBins)
    κBox   = zeros(T, EiBins, rhoBins, radBins)
    SBox   = zeros(T, EiBins, rhoBins, radBins)
    lnRoss = zeros(T, EiBins, rhoBins)

    B  = zeros(T, radBins)
    δB = zeros(T, radBins)
    #B_tot ::Float32 = 0.0

    bnu ::Float32  = 0.0
    dbnu ::Float32 = 0.0
    b::Int = 0

    ρ = exp.(eos.lnRho)
    Temp = exp.(eos.lnT)

    rtest = argmin(abs.(ρ .- test_r))
    ttest = argmin(abs.(eos.lnT[:, rtest] .- test_t))

    # None of the bins is allowed to be empty as this stage
    @assert all(binning .!= 0)

    @inbounds for j in eachindex(eos.lnRho)
        @inbounds for i in eachindex(eos.lnEi)
            B  .= 0.0
            δB .= 0.0
            #B_tot = 0.0
            @inbounds for k in eachindex(opacities.λ)
                b    = binning[k]
                bnu  = Bν( opacities.λ[k], Temp[i, j])
                dbnu = δBν(opacities.λ[k], Temp[i, j]) 

                χBox[ i, j, b]  += weights[k] .* opacities.κ[i, j, k]                            .* bnu  #Bν( opacities.λ[k], T[i, j])   
                χRBox[i, j, b]  += weights[k] .* 1.0 ./ opacities.κ[i, j, k]                     .* dbnu #δBν(opacities.λ[k], T[i, j]) 
                κBox[ i, j, b]  += weights[k] .* (opacities.κ[i, j, k] .- scattering.κ[i, j, k]) .* bnu  #Bν( opacities.λ[k], T[i, j])   
                SBox[ i, j, b]  += weights[k] .*                                                    bnu  #Bν( opacities.λ[k], T[i, j])  
                B[ b]           += weights[k] .*                                                    bnu  #Bν( opacities.λ[k], T[i, j]) 
                δB[b]           += weights[k] .*                                                    dbnu #δBν(opacities.λ[k], T[i, j]) 
            end
            #B_tot = sum(δB)
            #lnRoss[i, j]     = log( 1.0 / (sum(χRBox, dims=3) / B_tot))
            χBox[i, j, :]  ./= B
            κBox[i, j, :]  ./= B
            χRBox[i, j, :] ./= δB
        end

        if j == rtest
            @info "TestBin: $(test_bin), rho: $(ρ[j])"
            @info "χ: $(log(χBox[ttest, j, test_bin])), χR: $(log(1. /χRBox[ttest, j, test_bin])), simple mean: $(log(mean(χBox[ttest, j, test_bin])))"
            @info "wthin: $(exp(T(-1.5e7) .* T(1.0) ./χRBox[ttest, j, test_bin] *ρ[j]))"
        end

        χBox[:,  j, :] .*= ρ[j]
        κBox[:,  j, :] .*= ρ[j]
        χRBox[:, j, :] ./= ρ[j]
    end

    wthin  = exp.(T(-1.5e7) .* T(1.0) ./χRBox)
    wthick = T(1.0) .- wthin

    # Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    # Opacity
    opacity_table = remove_from_thin ? 
                    log.(wthin .* κBox .+ wthick .* (T(1.0) ./ χRBox)) : 
                    log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox))

    # ϵ table
    ϵ_table = log.(κBox ./ χBox)

    # Source function table --> S                = x/x+s * B + s/x+s * J 
    #                       --> thermal emission = x/x+s * B
    S_table = log.(SBox)
    #S_table = log.(κBox ./ χBox .* SBox)

    RegularOpacityTable(opacity_table, ϵ_table, S_table, collect(T, 1:radBins), false)
end


function box_integrated_v3(binning, weights, eos, opacities; remove_from_thin=false)
    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    EiBins  = length(eos.lnEi)
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, EiBins, rhoBins, radBins)    
    χRBox  = zeros(T, EiBins, rhoBins, radBins)
    κBox   = zeros(T, EiBins, rhoBins, radBins)
    SBox   = zeros(T, EiBins, rhoBins, radBins)

    B  = zeros(T, radBins)
    δB = zeros(T, radBins)

    bnu  ::Float32  = 0.0
    dbnu ::Float32  = 0.0
    b    ::Int      = 0

    ρ    = exp.(eos.lnRho)
    Temp = exp.(eos.lnT)

    rtest = argmin(abs.(ρ .- test_r))
    ttest = argmin(abs.(eos.lnT[:, rtest] .- test_t))

    # None of the bins is allowed to be empty as this stage
    @assert all(binning .!= 0)

    @inbounds for j in eachindex(eos.lnRho)
        @inbounds for i in eachindex(eos.lnEi)
            B  .= 0.0
            δB .= 0.0
            @inbounds for k in eachindex(opacities.λ)
                b    = binning[k]
                bnu  = Bν( opacities.λ[k], Temp[i, j])
                dbnu = δBν(opacities.λ[k], Temp[i, j]) 

                χBox[ i, j, b]  += weights[k] .* opacities.κ[i, j, k]                            .* bnu  #Bν( opacities.λ[k], T[i, j])   
                χRBox[i, j, b]  += weights[k] .* 1.0 ./ opacities.κ[i, j, k]                     .* dbnu #δBν(opacities.λ[k], T[i, j]) 
                SBox[ i, j, b]  += weights[k] .*                                                    bnu  #Bν( opacities.λ[k], T[i, j])  
                B[ b]           += weights[k] .*                                                    bnu  #Bν( opacities.λ[k], T[i, j]) 
                δB[b]           += weights[k] .*                                                    dbnu #δBν(opacities.λ[k], T[i, j]) 
            end
            χBox[i, j, :]  ./= B
            χRBox[i, j, :] ./= δB
        end

        if j == rtest
            @info "TestBin: $(test_bin), rho: $(ρ[j])"
            @info "χ: $(log(χBox[ttest, j, test_bin])), χR: $(log(1. /χRBox[ttest, j, test_bin])), simple mean: $(log(mean(χBox[ttest, j, test_bin])))"
            @info "wthin: $(exp(T(-1.5e7) .* T(1.0) ./χRBox[ttest, j, test_bin] *ρ[j]))"
        end

        χBox[:,  j, :] .*= ρ[j]
        χRBox[:, j, :] ./= ρ[j]
    end

    κBox  .= χBox
    wthin  = exp.(T(-1.5e7) .* T(1.0) ./χRBox)
    wthick = T(1.0) .- wthin

    ## Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    ## Opacity
    opacity_table = remove_from_thin ? 
                    log.(wthin .* κBox .+ wthick .* (T(1.0) ./ χRBox)) : 
                    log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox))

    ## ϵ table
    ϵ_table = log.(κBox ./ χBox)

    ## Source function table --> S                = x/x+s * B + s/x+s * J 
    ##                       --> thermal emission = x/x+s * B
    S_table = log.(SBox)
    #S_table = log.(κBox ./ χBox .* SBox)

    RegularOpacityTable(opacity_table, ϵ_table, S_table, collect(T, 1:radBins), false)
end

function clear!(args...)
    for arg in args
        arg .= 0.0
    end;
end

clear!(b::BinningStorage) = begin
    clear!([getfield(b, f) for f in fieldnames(typeof(b))]...)
    @assert all(store.wT .== 0.0)
end


"""
    tabulate(binning, weights, eos, opacities; kwargs...)

Tabulate opacity and source function for the given binning. Calls the currently
recommended version of box_integrated, which is box_integrated_v3, which is 
the fastest implementation, however without scattering at the moment to save memory.
"""
tabulate(args...; kwargs...) = box_integrated_v3(args...; kwargs...)

###########################################################################

## Debugging things
const test_bin = 4
const test_t   = 4843.07
const test_r   = 2.38802e-8