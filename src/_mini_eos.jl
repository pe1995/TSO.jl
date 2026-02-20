# ============================================================================
# MiniOpacityTable -- Opacity table without source function
# ============================================================================

"""
    MiniOpacityTable(opacity::RegularOpacityTable)

Create a reduced opacity table. This effectively does not store the source function, 
which saves disk space. However, this means that any lookup involving the source function
will be replaced by a simple call to the Planck function.
"""
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
reload(s_loc::Type{MiniOpacityTable}, path::String; binned) = begin
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
# Modified lookup function
# ============================================================================

function lookup(eos::E, s::MiniOpacityTable, v::Symbol, lnρ::AbstractFloat, lnT::AbstractFloat, i::Int) where {E1<:EoSTable, E2<:AxedEoS, E<:Union{E1, E2}}
    if v == :src
        return Bλ_fast(s.opacity.λ[i], exp(lnT)) * s.weights[i]
    elseif v == :dS_dT
        return dBdTλ_fast(s.opacity.λ[i], exp(lnT)) * s.weights[i]
    elseif v == :κ
        return lookup(eos, s.opacity, :κ, lnρ, lnT, i) * exp(lnρ)
    else
        return lookup(eos,s.opacity, v, lnρ, lnT, i)
    end
end

"""
    lookup_variable(s::MiniOpacityTable, v::Symbol)

Return a lookup function for the given variable. The function will take the following arguments:
v==:src 
    - T: Temperature
    - i: Index of the wavelength bin
v==:dS_dT
    - T: Temperature
    - i: Index of the wavelength bin
v==:κ
    - eos: EoS table
    - lnρ: Logarithm of density
    - lnT: Logarithm of temperature
    - i: Index of the wavelength bin
v==:other
    - eos: EoS table
    - lnρ: Logarithm of density
    - lnT: Logarithm of temperature
    - i: Index of the wavelength bin
"""
lookup_variable(s::MiniOpacityTable, v::Symbol) = begin
    if v == :src
        return (T, i) -> Bλ_fast(s.opacity.λ[i], T) * s.weights[i]
    elseif v == :dS_dT
        return (T, i) -> dBdTλ_fast(s.opacity.λ[i], T) * s.weights[i]
    elseif v == :κ
        return (eos, lnρ, lnT, i) -> lookup(eos, s.opacity, :κ, lnρ, lnT, i) * exp(lnρ)
    else
        return (eos, lnρ, lnT, i) -> lookup(eos,s.opacity, v, lnρ, lnT, i)
    end
end

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