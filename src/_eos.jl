# ============================================================================
# EoS related physics
# ============================================================================

# ============================================================================
# Planck function
# ============================================================================

function Bλ(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0)
    B[T .<= 1e-1] .= 0.0

    B
end

Bλ(λ::AbstractVector, T::AbstractVector) = begin
    B = zeros(eltype(T), length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= Bν(λ[i], T)
    end

    B
end

Bλ(λ::AbstractVector, T::AbstractMatrix) = begin
    B = zeros(eltype(T), size(T)..., length(λ))
    for i in axes(B, 3)
        B[:, :, i] .= Bν(λ[i], T)
    end

    B
end

function δBλ(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 * hc_k * exp(hc_k / (Λ*T)) / Λ^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2
    B[T .<= 1e-1] .= 0.0

    B
end

δBλ(λ::AbstractVector, T::AbstractVector) = begin
    B = zeros(eltype(T), length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= δBν(λ[i], T)
    end

    B
end

δBλ(λ::AbstractVector, T::AbstractMatrix) = begin
    B = zeros(eltype(T), size(T)..., length(λ))
    for i in axes(B, 3)
        B[:, :, i] .= δBν(λ[i], T)
    end

    B
end

@inline function δBλ(λ::A, T::A) where {A<:AbstractFloat}
    Λ = λ * aa_to_cm
    v = twohc2 * hc_k * exp(hc_k / (Λ*T)) / Λ^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2

    if isnan(v) || v<1e-30
        A(1e-30)
    else
        v
    end
end

@inline function Bλ(λ::A, T::A) where {A<:AbstractFloat}
    Λ = λ * aa_to_cm
    v = twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0)

    if isnan(v) || v<1e-30
        A(1e-30)
    else
        v
    end
end

Bν(args...; kwargs...)  = Bλ(args...; kwargs...)
δBν(args...; kwargs...) = δBλ(args...; kwargs...)

@inline Bν!(B::Ref{A}, λ::A, T::A) where {A<:AbstractFloat} = begin
    B[] = twohc2 /(λ * aa_to_cm)^5 /(exp(hc_k / ((λ * aa_to_cm)*T)) - 1.0)

    B[] = if isnan(B[]) || B[]<1e-30
        A(1e-30)
    else
        B[]
    end

    B[]
end

@inline δBν!(B::Ref{A}, λ::A, T::A) where {A<:AbstractFloat} = begin
    B[] = twohc2 * hc_k * exp(hc_k / ((λ * aa_to_cm)*T)) / (λ * aa_to_cm)^6 / T^2 / (exp(hc_k / ((λ * aa_to_cm)*T))-1)^2

    B[] = if isnan(B[]) || B[]<1e-30
        A(1e-30)
    else
        B[]
    end

    B[]
end

# ============================================================================
# compute Rosseland opacity
# ============================================================================

"""
    rosseland_opacity(tab; weights=tab.opacity.weights)

Integrate the rosseland opacity from the monochromatic opacity of the table.
"""
rosseland_opacity(tab::EquationOfState; kwargs...) = begin
    lnRoss = similar(tab[:lnRoss])
    rosseland_opacity!(lnRoss, tab; kwargs...)
    lnRoss
end

_rosseland_opacityX!(lnRoss, lnRho, axis_val, l, weights, opacity, T) = begin
    T_op = eltype(opacity)

    invT = zeros(T_op, size(T))
    @inbounds for i in eachindex(T)
        invT[i] = T[i] > 1e-1 ? 1.0 / T[i] : 0.0
    end

    nl = length(l)
    nt = length(axis_val)
    nr = length(lnRho)

    wk_c1c2_arr = zeros(T_op, nl)
    c2_arr = zeros(T_op, nl)
    @inbounds for k in 1:nl
        wk = T_op(weights[k])
        lk = T_op(l[k])
        Λ  = lk * T_op(aa_to_cm)
        c2 = T_op(hc_k) / Λ
        c1 = T_op(twohc2) / (Λ^5)

        c2_arr[k] = c2
        wk_c1c2_arr[k] = wk * c1 * c2
    end

    n_chunks = Threads.nthreads()
    chunk_size = cld(nl, n_chunks)
    chunks = Iterators.partition(1:nl, chunk_size) |> collect

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            l_ross = zeros(T_op, nt, nr)
            l_b    = zeros(T_op, nt, nr)

            @inbounds for k in chunk
                wk_c1c2 = wk_c1c2_arr[k]
                c2      = c2_arr[k]

                for j in 1:nr
                    @turbo for i in 1:nt
                        iT = invT[i, j]

                        x = c2 * iT
                        E = exp(x)

                        E_minus_1 = ifelse(iT > 0.0, E - 1.0, 1.0)

                        dbnuChunk_wk = (wk_c1c2 * iT * iT * E) / (E_minus_1 * E_minus_1)

                        val = ifelse((dbnuChunk_wk == dbnuChunk_wk) & (dbnuChunk_wk >= 1e-30), dbnuChunk_wk, T_op(1e-30))
                        val = ifelse(iT > 0.0, val, T_op(1e-30))

                        op = opacity[i, j, k]

                        val_valid = (op > 0.0)
                        vop = ifelse(val_valid, val / op, zero(T_op))
                        vb  = ifelse(val_valid, val, zero(T_op))

                        l_ross[i, j] += vop
                        l_b[i, j]    += vb
                    end
                end
            end
            return (l_ross, l_b)
        end
    end

    results = fetch.(tasks)

    ross_sum = zeros(T_op, nt, nr)
    b_sum    = zeros(T_op, nt, nr)

    for (r, b) in results
        ross_sum .+= r
        b_sum    .+= b
    end

    lnRoss .= 0.0
    @inbounds for j in 1:nr
        for i in 1:nt
            b = b_sum[i, j]
            if b > 0.0
                r = ross_sum[i, j] / b
                lnRoss[i, j] = (r == 0) ? T_op(NaN) : -log(r)
            else
                lnRoss[i, j] = T_op(NaN)
            end
        end
    end

    return lnRoss
end

"""
    rosseland_opacity!(lnRoss, tab; weights=tab.opacity.weights)

Integrate the rosseland opacity from the monochromatic opacity of the table.
"""
function rosseland_opacity!(lnRoss, tab::EquationOfState; weights=tab.opacity.weights)
    T = zeros(eltype(tab.lnT), size(tab)...)
    if ndims(tab.lnT) == 2
        T .= exp.(tab.lnT)
    else
        puff_up!(T, exp.(tab.lnT))
    end

    @optionalTiming rosseland_time _rosseland_opacityX!(
        lnRoss, tab.lnRho, tab.energy_axis.values, wavelength(tab), weights, tab[:κ], T
    )
end


# ============================================================================
# Optical depth related functions
# ============================================================================

_optical_depth!(ρκ, τ, z, λ, lnρ, lnE, eos, binned) = begin
    cRho, cE = weights(eos, lnρ, lnE; extrapolate=true)
    @inbounds for i in eachindex(λ)
        sample!((ρκ,), eos, (:κ,), cRho, cE, i)
        ρκ .= binned ? ρκ : exp.(lnρ) .* ρκ

        @inbounds for j in eachindex(z)
            if j==1
                τ[1, i] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
            else
                τ[j, i] = τ[j-1, i] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
            end
        end
    end

    τ
end

_optical_depthX!(τ, z, λ, lnρ, lnE, eos, binned) = begin
    chunks = Iterators.partition(eachindex(λ), cld(length(λ), Threads.nthreads())) |> collect

    tasks = map(chunks) do chunk
        Threads.@spawn _optical_depth_chunk!(chunk, z, λ, lnρ, lnE, eos, binned)
    end
    results = fetch.(tasks)

    for (i, chunk) in enumerate(chunks)
        τ[:, chunk] .= results[i]
    end
end

_optical_depth_chunk!(chunk, z, λ, lnρ, lnE, eos, binned) = _optical_depth_chunk!(
    zeros(eltype(eos[:κ]), length(z), length(chunk)), zeros(eltype(eos[:κ]), length(z)),
    chunk, z, λ, lnρ, lnE, eos, binned
)

_optical_depth_chunk!(τ, ρκ, chunk, z, λ, lnρ, lnE, eos, binned) = begin
    ρ = exp.(lnρ)
    cRho, cE = weights(eos, lnρ, lnE; extrapolate=true)

    @inbounds for (i, k) in enumerate(chunk)
        sample!((ρκ,), eos, (:κ,), cRho, cE, k)
        ρκ .= binned ? ρκ : ρ .* ρκ

        @inbounds for j in eachindex(z)
            if j==1
                τ[1, i] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
            else
                τ[j, i] = τ[j-1, i] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
            end
        end
    end

    τ
end

"""
    optical_depth(tab, model; binned=false)

Compute monochromatic and rosseland optical depth scales of the model
based on the opacity of the table.
"""
function optical_depth(eos::EquationOfState, model::AbstractModel; binned=false)
    model = TSO.flip(model, depth=true)

    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos[:κ])
    λ = wavelength(eos)

    τ_λ = zeros(T, length(lnρ), length(λ))
    τ_ross = zeros(T, length(lnρ))

    @optionalTiming optical_depth_time if Threads.nthreads() <= 1
        ρκ = zeros(T, length(lnρ))
        _optical_depth!(ρκ, τ_λ, z, λ, lnρ, lnE, eos, binned)
    else
        _optical_depthX!(τ_λ, z, λ, lnρ, lnE, eos, binned)
    end

    τ_ross .= rosseland_optical_depth(eos, model)

    τ_λ[τ_λ .> 1e30] .= 1e30
    τ_λ[τ_λ .< 1e-30] .= 1e-30

    τ_ross[τ_ross .> 1e30] .= 1e30
    τ_ross[τ_ross .< 1e-30] .= 1e-30

    return τ_ross, τ_λ
end

_rosseland_κ(eos, lnρ, lnE) = exp.(sample(eos, (:lnRoss,), lnρ, lnE; extrapolate=true)[1])

"""
    rosseland_optical_depth(tab, model)

Rosseland optical depth of the model, from the Rosseland opacity of the EoS.
"""
function rosseland_optical_depth(eos::EquationOfState, model::AbstractModel)
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos[:lnRoss])

    τ_ross = zeros(T, length(lnρ))
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    κ  .= _rosseland_κ(eos, lnρ, lnE)
    ρκ .= exp.(lnρ) .* κ
    for j in eachindex(z)
        if j==1
            τ_ross[1] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
        else
            τ_ross[j] = τ_ross[j-1] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
        end
    end

    return τ_ross
end

"""
    rosseland_depth(tab, model)

Geometrical depth of a model given on the Rosseland optical depth scale.
"""
function rosseland_depth(eos::EquationOfState, model::AbstractModel)
    τ, lnρ, lnE = model.τ, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos[:lnRoss])

    z      = zeros(T, length(lnρ))
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    κ  .= _rosseland_κ(eos, lnρ, lnE)
    ρκ .= exp.(lnρ) .* κ
    for j in eachindex(z)
        if j==1
            z[1] = 0 + (τ[2] - τ[1]) * 0.5 * (1.0/ρκ[j])
        else
            z[j] = z[j-1] + (τ[j] - τ[j-1]) * 0.5 * (1.0/ρκ[j] + 1.0/ρκ[j-1])
        end
    end

    return z
end

# ============================================================================
# formation height
# ============================================================================

"""
    formation_height(model, eos, τ_ross, τ_λ)

Compute the formation height + opacity, i.e. the rosseland optical depth
where the monochromatic optical depth is 1.
"""
function formation_height(model::AbstractModel, eos::EquationOfState, τ_ross, τ_λ)
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos[:κ])

    z   = Base.convert.(T, z)
    lnρ = Base.convert.(T, lnρ)
    lnE = Base.convert.(T, lnE)

    rosseland_depth = zeros(T, size(τ_λ, 2))
    opacity_depth   = zeros(T, size(τ_λ, 2))

    lRoss = log10.(τ_ross)
    lλ    = log10.(τ_λ)

    λ = wavelength(eos)

    r_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnρ, extrapolation_bc=Flat())
    T_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnE, extrapolation_bc=Flat())

    if Threads.nthreads() > 1
        @info "Computing formation height with $(Threads.nthreads()) threads."
    end

    @optionalTiming formation_height_time begin
        chunks = Iterators.partition(eachindex(λ), cld(length(λ), Threads.nthreads())) |> collect

        tasks = map(chunks) do chunk
            Threads.@spawn _formation_height_chunk(chunk, τ_λ, lλ, lRoss, r_ross, T_ross, eos)
        end
        results = fetch.(tasks)
        rossChunk = getindex.(results, 1)
        opaChunk = getindex.(results, 2)

        for (i, chunk) in enumerate(chunks)
            rosseland_depth[chunk] .= rossChunk[i]
            opacity_depth[chunk] .= opaChunk[i]
        end
    end

    exp10.(rosseland_depth), opacity_depth
end

_formation_height_chunk(chunk, τ_λ, lλ, lRoss, r_ross, T_ross, eos) = begin
    T = eltype(eos[:κ])
    rosseland_depth = zeros(T, length(chunk))
    opacity_depth = zeros(T, length(chunk))
    t_mono = zeros(T, size(τ_λ, 1))

    @inbounds for (i, k) in enumerate(chunk)
        t_mono .= @view(lλ[:, k])
        rosseland_depth[i] = linear_interpolation(
            Interpolations.deduplicate_knots!(t_mono, move_knots=true),
            lRoss,
            extrapolation_bc=Line()
        )(0.0)

        opacity_depth[i] = sample(eos, (:κ,), r_ross(rosseland_depth[i]), T_ross(rosseland_depth[i]), k; extrapolate=true)[1]
    end

    return [rosseland_depth, opacity_depth]
end

const rosseland_time = Ref(false)
const optical_depth_time = Ref(false)
const formation_height_time = Ref(false)

include("_interpolations.jl")

# ============================================================================
# Timer setup
# ============================================================================
const rosseland_time = Ref(false)
const optical_depth_time = Ref(false)
const formation_height_time = Ref(false)
