# ============================================================================
# MiniOpacityTable -- Opacity table without source function
# ============================================================================

struct MiniOpacityTable <: AbstractRegularTable
    opacity::RegularOpacityTable
    binned::Bool
    weights
end

get_wavelength(path::String) = begin
    fid = h5open(path, "r")
    λ = read(fid["λ"])
    close(fid)
    λ
end

"""
    reload(s::Type{MiniOpacityTable}, path::String)

Skip loading of the source function.
"""
reload(s_loc::Type{MiniOpacityTable}, path::String; binned=false) = begin
    fid   = HDF5.h5open(path, "r")

    λ = get_from_hdf5(Array, fid, :λ)
    κ = get_from_hdf5(Array, fid, :κ, mmap=true)
    κ_ross = get_from_hdf5(Array, fid, :κ_ross)
    optical_depth = get_from_hdf5(Bool, fid, :optical_depth)
    src = zeros(eltype(κ), 1, 1, 1)

    close(fid)

    o = RegularOpacityTable(κ, κ_ross, src, λ, optical_depth)
    MiniOpacityTable(o, binned, binned ? ones(length(λ)) : ω_midpoint(o))
end

# ============================================================================
# Accessor
# ============================================================================

opacity(opa::MiniOpacityTable) = opa.opacity
wavelength(opa::MiniOpacityTable) = wavelength(opa.opacity)

# ============================================================================
# Fast Planck function and its derivative
# ============================================================================

@inline function dBdTλ_fast(λ::AbstractFloat, T::AbstractFloat)
    T <= 1e-1 && return 0.0

    @fastmath begin
        Λ = λ * aa_to_cm

        invT = 1.0 / T

        c2 = hc_k / Λ
        x = c2 * invT

        E = exp(x)
        E_minus_1 = E - 1.0

        c1 = twohc2 / Λ^5
        return (c1 * x * invT * E) / (E_minus_1 * E_minus_1)
    end
end

@inline function Bλ_fast(λ::AbstractFloat, T::AbstractFloat)
    T <= 1e-1 && return 0.0

    @fastmath begin
        Λ = λ * aa_to_cm

        invT = 1.0 / T
        c2 = hc_k / Λ

        E = exp(c2 * invT)

        c1 = twohc2 / Λ^5

        return c1 / (E - 1.0)
    end
end
