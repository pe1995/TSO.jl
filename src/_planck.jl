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
