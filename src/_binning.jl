# ============================================================================
# Binned opacities
# ============================================================================

struct BinnedOpacities{O<:RegularOpacityTable, F<:AbstractFloat, Nϵ, Nw} <:AbstractBinnedOpacities
    opacities ::O
    ϵ         ::Array{F, Nϵ}
    wthin     ::Array{F, Nw}
end

# ============================================================================
# Wavelength integration weights
# ============================================================================

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

include("_advanced_binning.jl")

# ============================================================================
# Extensions of EoS functions for binned opacities
# ============================================================================

optical_depth(aos::E, opacities::B, model::AbstractModel) where {E<:AxedEoS, B<:BinnedOpacities} = optical_depth(aos, opacities.opacites, model, binned=true)

rosseland_optical_depth(eos::E,        opacities::BinnedOpacities, model::AbstractModel; binned=false) where {E<:AxedEoS} = rosseland_optical_depth(eos, opacities.opacities, model; binned=true)
rosseland_optical_depth(eos::EoSTable, opacities::BinnedOpacities, model::AbstractModel; binned=false) = rosseland_optical_depth(@axed(eos), opacities.opacities, model; binned=true)

formation_height(model::AbstractModel, eos::E, opacities::BinnedOpacities, τ_ross, τ_λ) where {E<:AxedEoS     } = formation_height(model, eos, opacities.opacities, τ_ross, τ_λ)
formation_height(model::AbstractModel, eos::E, opacities::BinnedOpacities, τ_ross, τ_λ) where {E<:OpacityTable} = formation_height(model, @axed(eos), opacities.opacities, τ_ross, τ_λ)

# ============================================================================
# Wrapper for normal functionality of binned tables
# ============================================================================

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
    _binned_from_regular!(opa_new, eos)
end

_binned_from_regular!(opa::SqOpacity, eos::SqEoS) = _binned_from_regular!(opa, @axed(eos))
_binned_from_regular!(opa::SqOpacity, eos::AxedEoS) = begin
    _, rr = meshgrid(eos)
    rr .= exp.(rr)
    weights = ω_midpoint(opa)

    Threads.@threads for i in eachindex(opa.λ)
        opa.κ[:, :, i] .*= rr
        opa.src[:, :, i] .*= weights[i]
    end

    BinnedOpacities(opa)
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

macro binned!(opa, eos)
    opa_l = esc(opa)
    eos_l = esc(eos)
    :(_binned_from_regular!($opa_l, $eos_l))
end

mean_intensity(aos::AxedEoS, opacity, model) = mean_intensity(aos, @binned(opacity, aos), model)
