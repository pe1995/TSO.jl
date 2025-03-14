#============================ Binning Interface ========================================#

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

struct BinnedOpacities{O<:RegularOpacityTable, F<:AbstractFloat, Nϵ, Nw} <:AbstractBinnedOpacities
    opacities ::O
    ϵ         ::Array{F, Nϵ}
    wthin     ::Array{F, Nw}
end

struct Beeck2012StaggerBins{F<:AbstractFloat} <:OpacityBins
    bin_edges ::Array{F,2}
end

"""
Stagger-style binning, but with horizontally stretched κ bins.
"""
struct SemiStaggerBins{F<:AbstractFloat} <:OpacityBins
    bin_edges ::Array{F,2}
end

"""
Bins based on the density distribution of points.
"""
struct DensityBins{F<:AbstractFloat} <:OpacityBins
    bin_edges ::Array{F,2}
end

"""
Kmeans Clustering bins.
"""
struct KmeansBins <:OpacityBins
    result
end

"""
Binning Quadrant with quadrant edges as arguments.
"""
struct Quadrant
    λ_lims
    κ_lims
    nbins
    stripes
end



EqualTabgenBins()   = EqualTabgenBins(Float64[])
UniformTabgenBins() = UniformTabgenBins(Float64[])



## Alias
TabgenBins = Union{<:EqualTabgenBins, <:UniformTabgenBins, <:CustomTabgenBins}
ClusterBins = Union{<:KmeansBins}
TabgenBinning(t::Type{<:OpacityBins}, args...; kwargs...)           = fill(t, args...; kwargs...)
StaggerBinning(t::Type{<:StaggerBins}, args...; kwargs...)          = fill(t, args...; kwargs...)
StaggerBinning(t::Type{<:Beeck2012StaggerBins}, args...; kwargs...) = fill(t, args...; kwargs...)
StaggerBinning(t::Type{<:SemiStaggerBins}, args...; kwargs...)      = fill(t, args...; kwargs...)
Co5boldBinning(t::Type{<:Co5boldBins}, args...; kwargs...)          = fill(t, args...; kwargs...)
DensityBinning(t::Type{<:DensityBins}, args...; kwargs...)          = fill(t, args...; kwargs...)
MURaMBinning(args...; kwargs...)                                    = fill(CustomTabgenBins, bin_edges=[-99.0,0.0,2.0,4.0,99.0], args...; kwargs...)
ClusterBinning(t::Type{<:ClusterBins}, args...; kwargs...)          = fill(t, args...; kwargs...)

assignments(a::Array{T, 1}) where {T} = a

Quadrant(; λ_lims, κ_lims, nbins, stripes=nothing) = Quadrant(λ_lims, κ_lims, nbins, stripes)
Quadrant(λ_lims, κ_lims, nbins; stripes=nothing) = Quadrant(λ_lims, κ_lims, nbins, stripes)



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


fill(::Type{Beeck2012StaggerBins}; kwargs...) = begin
    bin_edges = zeros(12, 4)

    bin_edges[1,  4] = -1.46; bin_edges[1,  3] = 9.00 ; bin_edges[1,  1] = 1.0   ; bin_edges[1,  2] = 3809.0
    bin_edges[2,  4] = -3.81; bin_edges[2,  3] = -1.46; bin_edges[2,  1] = 1.0   ; bin_edges[2,  2] = 3809.0
    bin_edges[3,  4] = -17.0; bin_edges[3,  3] =-3.81 ; bin_edges[3,  1] = 1.0   ; bin_edges[3,  2] = 3809.0
    bin_edges[4,  4] = -0.62; bin_edges[4,  3] = 9.0  ; bin_edges[4,  1] = 3809.0 ; bin_edges[4,  2] = 5624.
    bin_edges[5,  4] = -0.62; bin_edges[5,  3] = 9.0  ; bin_edges[5,  1] = 5624.0 ; bin_edges[5,  2] = 21612.
    bin_edges[6,  4] = -1.50; bin_edges[6,  3] = -0.62; bin_edges[6,  1] = 3809.0 ; bin_edges[6,  2] = 6426.
    bin_edges[7,  4] = -2.28; bin_edges[7,  3] = -1.50; bin_edges[7,  1] = 3809.0 ; bin_edges[7,  2] = 7109.
    bin_edges[8,  4] = -14.0; bin_edges[8,  3] = -2.28; bin_edges[8,  1] = 3809.0 ; bin_edges[8,  2] = 16465.
    bin_edges[9,  4] = -0.62; bin_edges[9,  3] = 9.0  ; bin_edges[9,  1] = 21612.0; bin_edges[9,  2] = 1000000.0
    bin_edges[10, 4] = -1.50; bin_edges[10, 3] = -0.62; bin_edges[10, 1] = 6426.0 ; bin_edges[10, 2] = 1000000.0
    bin_edges[11, 4] = -2.28; bin_edges[11, 3] = -1.5 ; bin_edges[11, 1] = 7109.0 ; bin_edges[11, 2] = 1000000.0
    bin_edges[12, 4] = -14.0; bin_edges[12, 3] = -2.28; bin_edges[12, 1] = 16465.0; bin_edges[12, 2] = 1000000.0

    bin_edges[:, 3] .*= -1
    bin_edges[:, 4] .*= -1
    bin_edges[:, 1] .= log10.(bin_edges[:, 1])
    bin_edges[:, 2] .= log10.(bin_edges[:, 2])
    Beeck2012StaggerBins(bin_edges)
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

function fill(::Type{<:SemiStaggerBins}; opacities, formation_opacity, Nbins=8, κ_bins=4, κ_edge=nothing)
    nλBins = Nbins - κ_bins
    logλ   = log10.(opacities.λ)

    formMin = minimum(formation_opacity)
    logλMin = minimum(logλ)

    λ_edges = zeros(nλBins+1)
    m  = isnothing(κ_edge) ? maximum(logλ) : maximum(logλ[formation_opacity .<= κ_edge])
    vm = isnothing(κ_edge) ? formation_opacity[argmax(logλ)] : κ_edge
    l  = minimum(logλ[formation_opacity .<= vm])
    r  = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    λ_edges .= collect(range(l, m, length=nλBins+1))

    ## reset the lower edge
    l  = logλMin
    m  = maximum(logλ)
    r  = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    λ_edges[1] = l
    λ_edges[end] = m

    ## Above the vm edge we collect everything in the remaining bins, splitting them equally
    k_edges = zeros(κ_bins+1)
    l = vm
    m = maximum(formation_opacity)
    #r = abs(m - l)
    #l = l - 0.01*r
    m = m + 0.01*r
    #l_bin_k_edge = vm

    # Split the part above the first bin equally
    k_splitted  = split_similar(sort(formation_opacity[formation_opacity.>l]), κ_bins)
    k_edges[1]  = l
    for b in 1:κ_bins
        k_edges[b+1] = last(k_splitted[b])
    end
    #k_edges .= collect(range(l, m, length=n_κ_bins+1))
    k_edges[end] = m

   
    bins = zeros(Nbins , 4)
    for i in 1:nλBins
        bins[i, 1] = λ_edges[i]
        bins[i, 2] = λ_edges[i+1]
        bins[i, 3] = formMin
        bins[i, 4] = vm
    end
    for i in 1:κ_bins
        bins[nλBins+i, 1] = λ_edges[1]
        bins[nλBins+i, 2] = λ_edges[end]
        bins[nλBins+i, 3] = k_edges[i]
        bins[nλBins+i, 4] = k_edges[i+1]
    end

   @show bins
   SemiStaggerBins(bins)
end

"""
Equal Tabgen bins based on opacities.
"""
function fill(::Type{<:EqualTabgenBins}; opacities, formation_opacity, Nbins=4, binsize=0.0, upper_limit=6.0, kwargs...)
    # Get the minimum and maximum opacity and decide on a step
    min_opacity  = minimum(formation_opacity)
    max_opacity  = isnothing(upper_limit) ? maximum(formation_opacity) : upper_limit

    r = abs(max_opacity - min_opacity)
    return if binsize==0.0
        if isnothing(upper_limit)
            EqualTabgenBins(collect(range(min_opacity-0.01*r, max_opacity+0.01*r, length=Nbins+1)))
        else
            b = collect(range(min_opacity-0.01*r, max_opacity, length=Nbins))
            append!(b, [maximum(formation_opacity) + 0.01*r])
            EqualTabgenBins(b)
        end
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
function fill(::Type{<:UniformTabgenBins}; opacities, formation_opacity, Nbins=4, line_bins=nothing, uniform_continuum=false, kwargs...)
    edges = if isnothing(line_bins)
        # sort the opacities and split the array in equal pieces
        sorted_opacities   = sort(formation_opacity)
        splitted_opacities = TSO.split_similar(sorted_opacities, Nbins)
        
        edges = [b[1] for b in splitted_opacities]
        append!(edges, splitted_opacities[end][end])

        edges
    else
        κ_edge = formation_opacity[argmax(opacities.λ)]

        cont_bins = Nbins - line_bins
        all_below = formation_opacity[formation_opacity .<= κ_edge]
        all_above = formation_opacity[formation_opacity .> κ_edge]

        sorted_opacities   = sort(all_below)
        edges = if !uniform_continuum
            splitted_opacities = TSO.split_similar(sorted_opacities, cont_bins)
            [b[1] for b in splitted_opacities]
        else
            collect(range(minimum(sorted_opacities), maximum(sorted_opacities), length=cont_bins+1))[1:end-1]
        end
    

        sorted_opacities   = sort(all_above)
        splitted_opacities = TSO.split_similar(sorted_opacities, line_bins)
        
        for b in splitted_opacities
            append!(edges, [b[1]])
        end

        append!(edges, splitted_opacities[end][end])

        edges
    end

    UniformTabgenBins(edges)
end

fill(::Type{<:ExactTabgenBins}; Nbins=4, dbox=0.5, kwargs...) = ExactTabgenBins(dbox, Nbins)

"""
stack 2 bins on top of each other, with the limit being 1 σ off the mean.
"""
function fill(::Type{DensityBins}; opacities, formation_opacity, λ_bins=4, κ_edge=nothing)
    nλBins = λ_bins
    logλ   = log10.(opacities.λ)

    formMin = minimum(formation_opacity)
    formMax = maximum(formation_opacity)
    rminmax = abs(formMax - formMin)
    formMax = formMax + 0.01*rminmax
    formMin = formMin - 0.01*rminmax

    logλMin = minimum(logλ)

    λ_edges = zeros(nλBins+1)
    m  = isnothing(κ_edge) ? maximum(logλ) : maximum(logλ[formation_opacity .<= κ_edge])
    vm = isnothing(κ_edge) ? formation_opacity[argmax(logλ)] : κ_edge
    l  = minimum(logλ[formation_opacity .<= vm])
    r  = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    λ_edges .= collect(range(l, m, length=nλBins+1))

    ## reset the lower edge
    l  = logλMin
    m  = maximum(logλ)
    r  = abs(m - l)
    l = l - 0.01*r
    m = m + 0.01*r
    λ_edges[1] = l
    λ_edges[end] = m

    λ_mask  = [(logλ .> λ_edges[i]) .& (logλ .< λ_edges[i+1]) for i in 1:nλBins]
    k_limit = [mean(formation_opacity[λ_mask[i]]) + 1.5*std(formation_opacity[λ_mask[i]]) for i in 1:nλBins]



    # Now we stack line bins on top
    bins = zeros(nλBins*2 , 4)
    for i in 1:nλBins
        bins[i, 1] = λ_edges[i]
        bins[i, 2] = λ_edges[i+1]
        bins[i, 3] = formMin
        bins[i, 4] = k_limit[i]

        bins[nλBins + i, 1] = λ_edges[i]
        bins[nλBins + i, 2] = λ_edges[i+1]
        bins[nλBins + i, 3] = k_limit[i]
        bins[nλBins + i, 4] = formMax
    end

    DensityBins(bins)
end

function fill(::Type{KmeansBins}; opacities, formation_opacity, Nbins=5, λ_split=nothing, κ_split=nothing, quadrants=nothing, stripes=false, kwargs...)
    if !isnothing(quadrants)
        masks = []
        clusters = []
        ll = log10.(opacities.λ)

        if sum([q.nbins for q in quadrants]) != Nbins
            @warn "The sum of all quadrants does not equal Nbins! Nbins will be $(sum([q.nbins for q in quadrants]))"
        end

        for q in quadrants
            nbins_split = q.nbins
            llims = q.λ_lims
            klims = q.κ_lims

            mask =  (formation_opacity .<= klims[2]) .& (formation_opacity .> klims[1]) .&
                    (ll .<= llims[2]) .& (ll .> llims[1])

            data = if isnothing(q.stripes) 
                hcat(ll[mask], formation_opacity[mask])'
            elseif q.stripes == :κ
                hcat(formation_opacity[mask])' 
            elseif q.stripes == :λ
                hcat(ll[mask])'
            else
                error("Please specify quadrant.stripes as `:κ`, `:λ`, or `nothing`")
            end

            cluster = kmeans(data, nbins_split; kwargs...)

            append!(clusters, [cluster])
            append!(masks, [mask])
        end
        
        as = zeros(Int, length(ll))
        for (i,mask) in enumerate(masks)
            max_as = maximum(as[.!mask], init=0.0)
            as[mask] = assignments(clusters[i]) .+ max_as
        end

        KmeansBins([as])
    elseif !isnothing(λ_split)
        nbins_split = λ_split[2]
        ll = log10.(opacities.λ)

        iend = findfirst(x->x>=λ_split[1], ll)
        data = stripes ? hcat(formation_opacity[iend+1:end])' : hcat(ll[iend+1:end], formation_opacity[iend+1:end])'
        clusters = kmeans(data, nbins_split; kwargs...)
    

        nbins_split = Nbins - λ_split[2]
        data = hcat(ll[1:iend], formation_opacity[1:iend])'
        clusters2 = kmeans(data, nbins_split; kwargs...)

        KmeansBins((clusters2, clusters))
    elseif !isnothing(κ_split)
        nbins_split = κ_split[2]
        ll = log10.(opacities.λ)

        mask = formation_opacity .<= κ_split[1]
        data = hcat(ll[mask], formation_opacity[mask])'
        clusters = kmeans(data, nbins_split; kwargs...)
    

        nbins_split = Nbins - κ_split[2]
        data = stripes ? hcat(formation_opacity[.!mask])' : hcat(ll[.!mask], formation_opacity[.!mask])'
        clusters2 = kmeans(data, nbins_split; kwargs...)

        as = zeros(Int, length(ll))
        as[mask]   = assignments(clusters)
        as[.!mask] = assignments(clusters2) .+ maximum(as[mask])

        KmeansBins([as])
    else
        data = stripes ? hcat(formation_opacity)' : hcat(log10.(opacities.λ), formation_opacity)'
        clusters = kmeans(data, Nbins; kwargs...)

        KmeansBins([clusters])
    end
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

function binning(b::B, opacities, formation_opacity, kwargs...) where {B<:Union{<:StaggerBins, <:Co5boldBins, <:Beeck2012StaggerBins, <:SemiStaggerBins}}
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

    @info "All wavelength points sorted in to bins? $(all(binning .!= 0))"
    binning
end

function binning(b::DensityBins, opacities, formation_opacity; splits=[], combine=[], kwargs...)
    bins_normal = binning(StaggerBins(b.bin_edges), opacities, formation_opacity, kwargs...)

    # Should bins be combined or splitted?
    combine_bins!(bins_normal, combine...)
    split_bins!(bins_normal, formation_opacity, splits...)

    #reset_bins(bins_normal)

    bins_normal
end

"""
Kmeans bin assignment
"""
binning(b::T, args...; kwargs...) where {T<:ClusterBins} = begin
    if length(b.result) == 1
        assignments(b.result |> first)
    else
        as = assignments.(b.result)
        for i in eachindex(as)
            i == 1 ? continue : nothing

            as[i] .+= maximum(as[i-1])
        end

        vcat(as...)
    end
end
    





"""
Reset the bins to a simple increasing structure.
"""
function reset_bins(binning)
    bin_assignment = []
    assigned_bins  = []

    for bin in binning
        if !(bin in assigned_bins)
            append!(assigned_bins, [bin])
            append!(bin_assignment, [length(bin_assignment)+1])
        end
    end

    bin_switch = Dict(a=>b for (a, b) in zip(assigned_bins, bin_assignment))

    binning_new = similar(binning)
    for (i, bin) in enumerate(binning)
        binning_new[i] = bin_switch[bin]
    end

    binning_new
end


function split_bins!(bins_normal, formation_opacity, splits...)
    nbins = maximum(unique(bins_normal))

    # Should bins be combined or splitted?
    for (i, split) in splits |> enumerate
        # bin to split
        i_bin_split = split[1]

        # in how many
        n_splits = split[2]
        
        # all points in that bin
        points_bin = bins_normal .== i_bin_split
        fopa_bin   = formation_opacity[points_bin]

        new_splits = split_similar(sort(fopa_bin), n_splits)
        k_edges = zeros(n_splits+1)
        k_edges[1] = minimum(fopa_bin)
        for b in 1:n_splits
            k_edges[b+1] = last(new_splits[b])
        end
        k_edges[end] = maximum(fopa_bin)

        l  = k_edges[1]
        m  = k_edges[end]
        r  = abs(m - l)
        #l = l - 0.01*r
        m = m + 0.01*r
        k_edges[1] = l
        k_edges[end] = m

        new_bin_numbers = [i_bin_split, (nbins + i for i in 1:n_splits-1)...]
        for i in eachindex(bins_normal)
            if !points_bin[i]
                continue
            end
            for j in 1:n_splits
                if (formation_opacity[i] >= k_edges[j]) & (formation_opacity[i] < k_edges[j+1])
                    bins_normal[i] = new_bin_numbers[j]
                end
            end
        end

        # make a new bin
        nbins = maximum(unique(bins_normal)) #nbins + (n_splits-1)
    end
end


function combine_bins!(bins_normal, combine...)
    for (i, comb) in combine |> enumerate
        new_bin_number = first(comb)

        for (i, i_bin) in bins_normal |> enumerate
            if i_bin in comb
                bins_normal[i] = new_bin_number
            end
        end 
    end
end


#======================= Wavelength integration weights ================================#

"""
Compute midpoint-weights from wavelenght array.
"""
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


"""
Compute midpoint-weights from wavelenght array between lower and upper boundary.
"""
function ω_midpoint(λ, lo, hi) 
    Δλ = zeros(eltype(λ), length(λ)+1)
    Δλ[2:end-1] = (λ[2:end] .- λ[1:end-1]) ./ 2
    Δλ[1]   = λ[1] - lo
    Δλ[end] = hi - λ[end]

    w  = similar(λ)
    for i in eachindex(λ)
        w[i] = Δλ[i] + Δλ[i+1]
    end

    w
end

ω_midpoint(opacities::OpacityTable) = ω_midpoint(to_cm(opacities.λ))





#====================== Computation of binned properties ===============================#

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

    BinnedOpacities(RegularOpacityTable(opacity_table, opacities.κ_ross, S_table, collect(T, 1:radBins), false), ϵ_table, wthin)
end

#=function box_integrated_v3(binning, weights, aos::E, opacities, scattering; remove_from_thin=false) where {E<:AxedEoS}
    eos = aos.eos
    eaxis = is_internal_energy(aos)
    
    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    AxBins  = length(aos.energy_axes.length)
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, AxBins, rhoBins, radBins)    
    χRBox  = zeros(T, AxBins, rhoBins, radBins)
    κBox   = zeros(T, AxBins, rhoBins, radBins)
    SBox   = zeros(T, AxBins, rhoBins, radBins)
    lnRoss = zeros(T, AxBins, rhoBins)
    Temp   = zeros(T, AxBins, rhoBins)
    if ndims(eos.lnT) == 2
        Temp .= exp.(eos.lnT)
    else
        puff_up!(Temp, exp.(eos.lnT))
    end

    B  = zeros(T, radBins)
    δB = zeros(T, radBins)
    #B_tot ::Float32 = 0.0

    bnu ::Float32  = 0.0
    dbnu ::Float32 = 0.0
    b::Int = 0

    ρ = exp.(eos.lnRho)

    #rtest = argmin(abs.(ρ .- test_r))
    #ttest = argmin(abs.(eos.lnT[:, rtest] .- test_t))

    # None of the bins is allowed to be empty as this stage
    @assert all(binning .!= 0)

    @inbounds for j in 1:rhoBins
        @inbounds for i in 1:AxBins
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

        #if j == rtest
        #    @info "TestBin: $(test_bin), rho: $(ρ[j])"
        #    @info "χ: $(log(χBox[ttest, j, test_bin])), χR: $(log(1. /χRBox[ttest, j, test_bin])), simple mean: $(log(mean(χBox[ttest, j, test_bin])))"
        #    @info "wthin: $(exp(T(-1.5e7) .* T(1.0) ./χRBox[ttest, j, test_bin] *ρ[j]))"
        #end

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

    BinnedOpacities(RegularOpacityTable(opacity_table, opacities.κ_ross, S_table, collect(T, 1:radBins), false), ϵ_table)
end=#

function box_integrated_v3(binning, weights, aos::E, opacities, scattering=nothing; transition_model=nothing) where {E<:AxedEoS}
    eos = aos.eos
    eaxis = is_internal_energy(aos)

    iscat = isnothing(scattering)
    remove_from_thin = !iscat

    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    AxBins  = aos.energy_axes.length
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, AxBins, rhoBins, radBins)    
    χRBox  = zeros(T, AxBins, rhoBins, radBins)
    χ_thin = similar(χRBox)
    κBox   = zeros(T, AxBins, rhoBins, radBins)
    SBox   = zeros(T, AxBins, rhoBins, radBins)
    Temp   = zeros(T, AxBins, rhoBins)
    if ndims(eos.lnT) == 2
        Temp .= exp.(eos.lnT)
    else
        puff_up!(Temp, exp.(eos.lnT))
    end

    B  = zeros(T, radBins)
    δB = zeros(T, radBins)

    bnu  ::Float32  = 0.0
    dbnu ::Float32  = 0.0
    b    ::Int      = 0

    ρ = exp.(eos.lnRho)

    ## The intensity can be computed from the transition model aswell, such that it is ready 
    ## for interpolation
    #J = isnothing(transition_model) ? nothing : @interpolated mean_intensity(aos, opacities, transition_model) :lnρ


    # None of the bins are allowed to be empty as this stage
    @assert all(binning .!= 0)

    @inbounds for j in 1:rhoBins
        @inbounds for i in 1:AxBins
            B  .= 0.0
            δB .= 0.0
            @inbounds for k in eachindex(opacities.λ)
                b    = binning[k]
                bnu  = Bν( opacities.λ[k], Temp[i, j]) #isnothing(J) ? Bν( opacities.λ[k], Temp[i, j]) : exp(evaluate_radiation(J, k, eos.lnRho[j]))
                dbnu = δBν(opacities.λ[k], Temp[i, j]) 

                χBox[ i, j, b]  += weights[k] .* opacities.κ[i, j, k]         .* bnu  
                χRBox[i, j, b]  += weights[k] .* 1.0 ./ opacities.κ[i, j, k]  .* dbnu 
                χ_thin[i,j, b]   = χRBox[i, j, b]
                SBox[ i, j, b]  += weights[k] .* bnu #Bν( opacities.λ[k], Temp[i, j]) #bnu    
                B[ b]           += weights[k] .* bnu  
                δB[b]           += weights[k] .* dbnu 
                κBox[ i, j, b]  += (!iscat) ? 
                                    weights[k] .* (opacities.κ[i, j, k] .- scattering.κ[i, j, k]) .* bnu : 
                                    weights[k] .* opacities.κ[i, j, k] .* bnu
            end
            χBox[i,   j, :] ./= B
            κBox[i,   j, :] ./= B
            χRBox[i,  j, :] ./= δB
            χ_thin[i, j, :] ./= δB
        end

        χBox[:,  j, :] .*= ρ[j]
        κBox[:,  j, :] .*= ρ[j]
        χRBox[:, j, :] ./= ρ[j]
    end

    χRBox[χRBox .<= 1e-30] .= 1e-30
    χRBox[χRBox .>= 1e30]  .= 1e30
    χBox[ χBox  .<= 1e-30] .= 1e-30
    χBox[ χBox  .>= 1e30]  .= 1e30
    SBox[SBox   .<= 1e-30] .= 1e-30
    SBox[SBox   .>= 1e30]  .= 1e30


    #@show minimum(χRBox) maximum(χRBox) minimum(χBox) maximum(χBox)

    κ_ross = T(1.0) ./χRBox

    wthin, wthick = if isnothing(transition_model)
        wthin_l = exp.(T(-1.5e7) .* κ_ross)
        wthick_l = T(1.0) .- wthin_l
        wthin_l, wthick_l
    else
        pure_rosseland = RegularOpacityTable(κ_ross, opacities.κ_ross, SBox, collect(T, 1:radBins), false)
        opacity_transition(aos, pure_rosseland, transition_model)
    end

    ## Use rosseland average where optically thick, Planck average where optically thin
    #return log.(wthin .* χBox .+ wthick .* (T(1.0) ./ χRBox)), log.(SBox), log.(κBox ./ χBox)
    #return log.(χBox), log.(SBox), log.(κBox ./ χBox)

    ## Opacity
    opacity_table = remove_from_thin ? 
                    wthin .* κBox .+ wthick .* κ_ross : 
                    wthin .* χBox .+ wthick .* κ_ross

    ## ϵ table
    ϵ_table = κBox ./ χBox

    ## Source function table --> S                = x/x+s * B + s/x+s * J 
    ##                       --> thermal emission = x/x+s * B
    S_table = SBox
    
    #S_table = log.(κBox ./ χBox .* SBox)
    #(BinnedOpacities(RegularOpacityTable(χBox, opacities.κ_ross, S_table, collect(T, 1:radBins), false), ϵ_table, wthin), 
    #BinnedOpacities(RegularOpacityTable(κ_ross, opacities.κ_ross, SBox, collect(T, 1:radBins), false)))

    BinnedOpacities(RegularOpacityTable(opacity_table, opacities.κ_ross, S_table, collect(T, 1:radBins), false), ϵ_table, wthin)
end

box_integrated_v3(binning, weights, eos::E, opacities, scattering; kwargs...) where {E<:RegularEoSTable} = box_integrated_v3(binning, weights, AxedEoS(eos), opacities, scattering; kwargs...)
box_integrated_v3(binning, weights, eos::E, opacities; kwargs...) where {E<:RegularEoSTable} = box_integrated_v3(binning, weights, AxedEoS(eos), opacities; kwargs...)







#= Faster version (v4) =#
#= 
This version uses precompilation and better ordering of loops
It specifically does not incluse SIMD, because it makes it slower for some reason 
=#

"""
    do_binningX!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
                        ρ, λ, binning, Temp, weights, κ)

Multithreaded version of do_binning. Uses partial chunk sums to avoid race conditions.
"""
function do_binningX!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
                        ρ, λ, binning, Temp, weights, κ, src, κ_scat; corr_χ=nothing, corr_S=nothing)
    rhoBins = size(χBox, 2)
    AxBins = size(χBox, 1)
    radBins = size(χBox, 3)

    # chucks for threading
    chunks = Iterators.partition(eachindex(λ), length(λ) ÷ Threads.nthreads()) |> collect

    begin
        tasks = map(chunks) do chunk
            Threads.@spawn _do_binning_chunks(chunk, rhoBins, AxBins, radBins, binning, weights, κ, κ_scat, src, Temp, λ; corr_χ=corr_χ, corr_S=corr_S)
        end
        results = fetch.(tasks)
        χBoxChunk = getindex.(results, 1)
        χRBoxChunk = getindex.(results, 2)
        χthinChunk = getindex.(results, 3)
        SBoxChunk = getindex.(results, 4)
        BChunk = getindex.(results, 5)
        δBChunk = getindex.(results, 6)
        κBoxChunk = getindex.(results, 7)
    end

    # sum over all chunks to get the total sum
    χBox .= sum(χBoxChunk)
    χRBox .= sum(χRBoxChunk)
    χ_thin .= sum(χthinChunk) 
    SBox .= sum(SBoxChunk)
    B .= sum(BChunk)
    δB .= sum(δBChunk)
    κBox .= sum(κBoxChunk)

    χBox   ./= B
    κBox   ./= B
    χRBox  ./= δB
    χ_thin ./= δB

    @inbounds for j in 1:rhoBins
        χBox[:,  j, :] .*= ρ[j]
        κBox[:,  j, :] .*= ρ[j]
        χRBox[:, j, :] ./= ρ[j]
    end

    χRBox[χRBox .<= 1e-30] .= 1e-30
    χRBox[χRBox .>= 1e30]  .= 1e30
    χBox[ χBox  .<= 1e-30] .= 1e-30
    χBox[ χBox  .>= 1e30]  .= 1e30
    κBox[ κBox  .<= 1e-30] .= 1e-30
    κBox[ κBox  .>= 1e30]  .= 1e30
    SBox[SBox   .<= 1e-30] .= 1e-30
    SBox[SBox   .>= 1e30]  .= 1e30
end

_do_binning_chunks(chunk, rhoBins, AxBins, radBins, binning, weights, κ, κ_scat, src, Temp, λ; corr_χ=nothing, corr_S=nothing) = begin
    bC = Ref{Int}(0)
    dbnC = Ref(eltype(κ)(0.0)) 
    bnC = Ref(eltype(κ)(0.0))
    χBoxChunk = zeros(eltype(κ), AxBins, rhoBins, radBins)
    χRBoxChunk = deepcopy(χBoxChunk)
    χthinChunk = deepcopy(χBoxChunk)
    SBoxChunk = deepcopy(χBoxChunk)
    BChunk = deepcopy(χBoxChunk)
    δBChunk = deepcopy(χBoxChunk)
    κBoxChunk = deepcopy(χBoxChunk)

    _do_binning_chunks_core!(
        binning, chunk, rhoBins, AxBins, weights, κ, κ_scat, src, Temp, λ, 
        bC, bnC, dbnC, 
        χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
        corr_χ, corr_S
    )

    return [
        χBoxChunk,
        χRBoxChunk,
        χthinChunk,
        SBoxChunk,
        BChunk,
        δBChunk,
        κBoxChunk
    ]
end






_do_binning_chunks_core!(
    binning, chunk, rhoBins, AxBins, weights, κ, κ_scat::Nothing, src, Temp, λ, 
    bC, bnC, dbnC, 
    χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
    corr_χ::Nothing, corr_S::Nothing
    ) = begin
    @inbounds for k in chunk
        bC[] = binning[k] 
        #(bC[] == 0) && continue
        @inbounds for j in 1:rhoBins
            @inbounds for i in 1:AxBins
                Bν!(bnC, λ[k], Temp[i, j]) 
                δBν!(dbnC, λ[k], Temp[i, j]) 

                χBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[]
                χRBoxChunk[i, j, bC[]] += weights[k] * 1.0 / κ[i, j, k] * dbnC[]
                χthinChunk[i, j, bC[]] = χBoxChunk[i, j, bC[]]
                SBoxChunk[i, j, bC[]]  += weights[k] * src[i, j, k] 
                BChunk[i, j, bC[]] += weights[k] * bnC[]
                δBChunk[i, j, bC[]] += weights[k] * dbnC[]
                κBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[]
            end
        end
    end
end

_do_binning_chunks_core!(
    binning, chunk, rhoBins, AxBins, weights, κ, κ_scat::Nothing, src, Temp, λ, 
    bC, bnC, dbnC, 
    χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
    corr_χ::T, corr_S::T
    ) where {T<:AbstractArray} = begin
    @inbounds for k in chunk
        bC[] = binning[k] 
        #(bC[] == 0) && continue
        @inbounds for j in 1:rhoBins
            @inbounds for i in 1:AxBins
                Bν!(bnC, λ[k], Temp[i, j]) 
                δBν!(dbnC, λ[k], Temp[i, j]) 

                χBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[] * corr_χ[i, j, k]
                χRBoxChunk[i, j, bC[]] += weights[k] * 1.0 / (κ[i, j, k]*corr_χ[i, j, k]) * dbnC[]
                χthinChunk[i, j, bC[]] = χBoxChunk[i, j, bC[]]
                SBoxChunk[i, j, bC[]]  += weights[k] * src[i, j, k] * corr_S[i, j, k]
                BChunk[i, j, bC[]] += weights[k] * bnC[]
                δBChunk[i, j, bC[]] += weights[k] * dbnC[]
                κBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[] * corr_χ[i, j, k]
            end
        end
    end
end

_do_binning_chunks_core!(
    binning, chunk, rhoBins, AxBins, weights, κ, κ_scat::O, src, Temp, λ, 
    bC, bnC, dbnC, 
    χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
    corr_χ::Nothing, corr_S::Nothing
    ) where {O<:AbstractArray} = begin
    @inbounds for k in chunk
        bC[] = binning[k] 
        #(bC[] == 0) && continue
        @inbounds for j in 1:rhoBins
            @inbounds for i in 1:AxBins
                Bν!(bnC, λ[k], Temp[i, j]) 
                δBν!(dbnC, λ[k], Temp[i, j]) 

                χBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[]
                χRBoxChunk[i, j, bC[]] += weights[k] * 1.0 / κ[i, j, k] * dbnC[]
                χthinChunk[i, j, bC[]] = χBoxChunk[i, j, bC[]]
                SBoxChunk[i, j, bC[]]  += weights[k] * src[i, j, k] 
                BChunk[i, j, bC[]] += weights[k] * bnC[]
                δBChunk[i, j, bC[]] += weights[k] * dbnC[]
                κBoxChunk[i, j, bC[]] += weights[k] * (κ[i, j, k] - κ_scat[i, j, k]) * bnC[]
            end
        end
    end
end

_do_binning_chunks_core!(
    binning, chunk, rhoBins, AxBins, weights, κ, κ_scat::O, src, Temp, λ, 
    bC, bnC, dbnC, 
    χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
    corr_χ::T, corr_S::T
    ) where {T<:AbstractArray, O<:AbstractArray} = begin
    @inbounds for k in chunk
        bC[] = binning[k] 
        #(bC[] == 0) && continue
        @inbounds for j in 1:rhoBins
            @inbounds for i in 1:AxBins
                Bν!(bnC, λ[k], Temp[i, j]) 
                δBν!(dbnC, λ[k], Temp[i, j]) 

                χBoxChunk[i, j, bC[]] += weights[k] * κ[i, j, k] * bnC[] * corr_χ[i, j, k]
                χRBoxChunk[i, j, bC[]] += weights[k] * 1.0 / (κ[i, j, k]*corr_χ[i, j, k]) * dbnC[]
                χthinChunk[i, j, bC[]] = χBoxChunk[i, j, bC[]]
                SBoxChunk[i, j, bC[]]  += weights[k] * src[i, j, k] * corr_S[i, j, k]
                BChunk[i, j, bC[]] += weights[k] * bnC[]
                δBChunk[i, j, bC[]] += weights[k] * dbnC[]
                κBoxChunk[i, j, bC[]] += weights[k] * (κ[i, j, k] - κ_scat[i, j, k]) * bnC[] * corr_χ[i, j, k]
            end
        end
    end
end




function box_integrated_v4(binning, weights, aos::E, opacities, scattering=nothing; 
    transition_model=nothing, corr_χ=nothing, corr_S=nothing) where {E<:AxedEoS}
    eos = aos.eos
    eaxis = is_internal_energy(aos)

    iscat = isnothing(scattering)
    remove_from_thin = !iscat

    if !remove_from_thin
        @warn "No scattering opacity provided!"
    end

    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    AxBins  = aos.energy_axes.length
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, AxBins, rhoBins, radBins)    
    χRBox  = zeros(T, AxBins, rhoBins, radBins)
    χ_thin = similar(χRBox)
    κBox   = zeros(T, AxBins, rhoBins, radBins)
    SBox   = zeros(T, AxBins, rhoBins, radBins)
    Temp   = zeros(T, AxBins, rhoBins)

    if ndims(eos.lnT) == 2
        Temp .= exp.(eos.lnT)
    else
        puff_up!(Temp, exp.(eos.lnT))
    end

    B  = zeros(T, AxBins, rhoBins, radBins)
    δB = zeros(T, AxBins, rhoBins, radBins)

    bnu  ::Float32  = 0.0
    dbnu ::Float32  = 0.0
    b    ::Int      = 0

    ρ = exp.(eos.lnRho)

    # None of the bins are allowed to be empty as this 
    if !all(binning .!= 0)
        @warn "N0 = $(count(binning .== 0)), $(count(binning .== 0)/length(binning)*100.0)"
        @assert all(binning .!= 0)
    end

    if Threads.nthreads() > 1
        @info "Binning opacities with $(Threads.nthreads()) threads."
    end

    # Do stuff
    scat = remove_from_thin ? scattering.κ  : nothing
    do_binningX!(
        B, δB, SBox, κBox, χBox, χRBox, χ_thin, ρ,
        opacities.λ, binning, Temp, weights, 
        opacities.κ, opacities.src, scat;
        corr_χ=corr_χ, corr_S=corr_S
    )

    κ_ross = T(1.0) ./χRBox

    wthin, wthick = if isnothing(transition_model)
        @warn "No transition model was provided to the binning. "*
              "The binnig will only be valid for the sun!"
        wthin_l = exp.(T(-1.5e7) .* κ_ross)
        wthick_l = T(1.0) .- wthin_l
        wthin_l, wthick_l
    else
        pure_rosseland = RegularOpacityTable(
            κ_ross, opacities.κ_ross, SBox, collect(T, 1:radBins), false
        )
        opacity_transition(aos, pure_rosseland, transition_model)
    end

    ## Opacity
    #opacity_table = remove_from_thin ? 
    #    wthin .* κBox .+ wthick .* κ_ross : 
    #    wthin .* χBox .+ wthick .* κ_ross
    #
    # if remove_from_thin => k = chi_tot - chi_scat, else k = chi_tot
    opacity_table = wthin .* κBox .+ wthick .* κ_ross 

    ## ϵ table
    ϵ_table = κBox ./ χBox

    ## Source function table --> S                = x/x+s * B + s/x+s * J 
    ##                       --> thermal emission = x/x+s * B
    S_table =  SBox #opacity_table ./ (wthin .* χBox .+ wthick .* κ_ross) .* SBox

    BinnedOpacities(
        SqOpacity(
            opacity_table, opacities.κ_ross, S_table,
            collect(T, 1:radBins), false
        ), 
        ϵ_table, 
        wthin
    )
end

box_integrated_v4(binning, weights, eos::E, opacities, scattering; kwargs...) where {E<:RegularEoSTable} = box_integrated_v4(binning, weights, AxedEoS(eos), opacities, scattering; kwargs...)
box_integrated_v4(binning, weights, eos::E, opacities; kwargs...) where {E<:RegularEoSTable} = box_integrated_v4(binning, weights, AxedEoS(eos), opacities; kwargs...)



#================= Extensions of EoS functions for binned opacities ====================#

"""
    opacity_transition(aos, opacities, model)

Compute the transition limit between free streaming and diffusion based on a model and binned opacities.
At this point one could either
    A) Compute the optical depth for the model and interpolate in density to the table values, might be weird
    because the table is much larger in density than the model. So there is extrapolation involved then.
    B) Use the approximation from Ludwig 1999, i.e.
        τ ≈ κross * Pg / g
For simplicity we choose B) for the time being.
"""
function opacity_transition(aos::A, opacities::O, model::AbstractModel) where {A<:AxedEoS, O<:OpacityTable}
    ## Compute the optical depth from the pure rosseland opacities
    #τ,_ = optical_depth(aos, @binned(opacities), model) --> then interpolate/extrapolate to table rho
    
    T = eltype(opacities.κ)

    ## Approximation from Ludwig 1999
    τ = zeros(T, size(opacities.κ)...)
    for k in eachindex(opacities.λ)
        for j in eachindex(DensityAxis(aos).values)
            τ[:, j, k] .= opacities.κ[:, j, k] ./ exp(aos.eos.lnRho[j]) .* exp.(aos.eos.lnPg[:, j]) ./ exp10(model.logg)
        end
    end

    ## τ contains the optical depth computed from the rosseland mean in each bin
    wthin = exp.(T(-2.0) .* τ)

    (wthin, T(1.0) .- wthin)
end

optical_depth(aos::E, opacities::B, model::AbstractModel) where {E<:AxedEoS, B<:BinnedOpacities} = optical_depth(aos, opacities.opacites, model, binned=true)

rosseland_optical_depth(eos::E,        opacities::BinnedOpacities, model::AbstractModel; binned=false) where {E<:AxedEoS} = rosseland_optical_depth(eos, opacities.opacities, model; binned=true)
rosseland_optical_depth(eos::EoSTable, opacities::BinnedOpacities, model::AbstractModel; binned=false) = rosseland_optical_depth(@axed(eos), opacities.opacities, model; binned=true)

formation_height(model::AbstractModel, eos::E, opacities::BinnedOpacities, τ_ross, τ_λ) where {E<:AxedEoS     } = formation_height(model, eos, opacities.opacities, τ_ross, τ_λ)
formation_height(model::AbstractModel, eos::E, opacities::BinnedOpacities, τ_ross, τ_λ) where {E<:OpacityTable} = formation_height(model, @axed(eos), opacities.opacities, τ_ross, τ_λ)




## Binned storage extensions (deprecated) 

function clear!(args...)
    for arg in args
        arg .= 0.0
    end;
end

clear!(b::BinningStorage) = begin
    clear!([getfield(b, f) for f in fieldnames(typeof(b))]...)
    @assert all(store.wT .== 0.0)
end





#======================== Convenience binning interface =================================#

"""
    tabulate(binning, weights, eos, opacities; kwargs...)

Tabulate opacity and source function for the given binning. Calls the currently
recommended version of box_integrated, which is box_integrated_v3, which is 
the fastest implementation, however without scattering at the moment to save memory.
"""
tabulate(args...; kwargs...) = @optionalTiming binning_time box_integrated_v4(args...; kwargs...)





#=============== Wrapper for normal functionality of binned tables =====================#

for_dispatch(eos::EoSTable, opacities::BinnedOpacities, args...) = for_dispatch(eos, 
                                                                        opacities.opacities.κ, 
                                                                        opacities.opacities.src, 
                                                                        opacities.ϵ, args...)

switch_energy(aos::A, opacities::BinnedOpacities, args...; kwargs...) where {A<:AxedEoS} = begin
    eos_e, opa_e = switch_energy(aos, opacities.opacities, args...; kwargs...)

    binned_e     = RegularOpacityTable(opacities.ϵ, opacities.opacities.κ_ross, opacities.opacities.src, opacities.opacities.λ, opacities.opacities.optical_depth)
    eos_e2, ϵ_e   = switch_energy(aos, binned_e, args...; kwargs...)

    eos_e, BinnedOpacities(opa_e, ϵ_e.κ)
end

complement(aos_old::A1, aos_new::A2, opacities::BinnedOpacities, args...; kwargs...) where {A1<:AxedEoS, A2<:AxedEoS} = begin
    opa_e = complement(aos_old, aos_new, opacities.opacities, args...; kwargs...)

    binned_e = RegularOpacityTable(opacities.ϵ, opacities.opacities.κ_ross, opacities.opacities.src, opacities.opacities.λ, opacities.opacities.optical_depth)
    ϵ_e  = complement(aos_old, aos_new, binned_e, args...; kwargs...)

    BinnedOpacities(opa_e, ϵ_e.κ)
end

uniform(aos::A, opacities::BinnedOpacities, args...; kwargs...) where {A<:AxedEoS} = begin
    eos_e, opa_e = complement(aos, opacities.opacities, args...; kwargs...)

    binned_e     = RegularOpacityTable(opacities.ϵ, opacities.opacities.κ_ross, opacities.opacities.src, opacities.opacities.λ, opacities.opacities.optical_depth)
    eos_e2, ϵ_e  = complement(aos, binned_e, args...; kwargs...)

    eos_e, BinnedOpacities(opa_e, ϵ_e.κ)
end

BinnedOpacities(opa::SqOpacity) = BinnedOpacities(opa, ones(eltype(opa.κ), size(opa.κ)...))
BinnedOpacities(opa::SqOpacity, ϵ::AbstractArray) = BinnedOpacities(opa, ϵ, ones(eltype(opa.κ), size(opa.κ)...))

_binned_from_regular(opa::SqOpacity) = BinnedOpacities(opa)
_binned_from_regular(opa::SqOpacity, eos::SqEoS)   = _binned_from_regular(opa, @axed(eos))
_binned_from_regular(opa::SqOpacity, eos::AxedEoS) = begin
    opa_new = deepcopy(opa)
    _, rr   = meshgrid(eos)
    rr     .= exp.(rr)
    for i in eachindex(opa.λ)
        opa_new.κ[:, :, i] .*= rr
    end
    
    BinnedOpacities(opa_new)
end


macro binned(opa)
    opa_l = esc(opa)
    :(_binned_from_regular($opa_l))
end

macro binned(opa, eos)
    opa_l = esc(opa)
    eos_l = esc(eos)
    :(_binned_from_regular($opa_l, $eos_l))
end

"""
Mean intensity from unbinned opacity table.
"""
mean_intensity(aos::AxedEoS, opacity, model) = mean_intensity(aos, @binned(opacity, aos), model)




#============================ Debugging things =========================================#

## Debugging things
const test_bin = 4
const test_t   = 4843.07
const test_r   = 2.38802e-8

const binning_time = Ref(false)

#=======================================================================================#