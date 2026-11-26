# ============================================================================
# Binned opacities
# ============================================================================

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

# ============================================================================
# Binning
# ============================================================================

const USE_BINNING_THREADS = Ref(true)

"""
    binning(W_assign, weights, tab; scattering=tab.scattering, transition_model, logg, corr_χ, corr_S)

Bin the monochromatic opacities of `tab` with the bin assignment `W_assign`.
Returns the binned table, which shares the EoS of `tab`. Its opacity holds the binned
opacity `κ`, source function `src`, photon destruction probability `ϵ` and the
optically thin weight `wthin`.
"""
function binning(W_assign::AbstractMatrix, weights, tab::EquationOfState; scattering=tab.scattering,
                 transition_model=nothing, logg=nothing, corr_χ=nothing, corr_S=nothing)

    radBins = size(W_assign, 2)

    remove_from_thin = !isnothing(scattering)

    if !remove_from_thin
        @warn "No scattering opacity provided!"
    end

    lnRho = tab.lnRho
    rhoBins = length(lnRho)
    AxBins  = tab.energy_axis.length
    T_type  = eltype(lnRho)

    ϵ_table       = zeros(T_type, AxBins, rhoBins, radBins)
    opacity_table = zeros(T_type, AxBins, rhoBins, radBins)
    S_table       = zeros(T_type, AxBins, rhoBins, radBins)

    χBox   = zeros(T_type, AxBins, rhoBins, radBins)
    χRBox  = zeros(T_type, AxBins, rhoBins, radBins)
    κ_ross = zeros(T_type, AxBins, rhoBins, radBins)
    χ_thin = zeros(T_type, AxBins, rhoBins, radBins)
    κBox   = zeros(T_type, AxBins, rhoBins, radBins)
    SBox   = zeros(T_type, AxBins, rhoBins, radBins)
    B      = zeros(T_type, AxBins, rhoBins, radBins)
    δB     = zeros(T_type, AxBins, rhoBins, radBins)

    Temp = zeros(T_type, AxBins, rhoBins)
    if ndims(tab.lnT) == 2
        Temp .= exp.(tab.lnT)
    else
        puff_up!(Temp, exp.(tab.lnT))
    end

    ρ = exp.(lnRho)
    scat = remove_from_thin ? scattering[:κ] : nothing

    _binning_X!(
        B, δB, SBox, κBox, χBox, χRBox, χ_thin, ρ,
        wavelength(tab), W_assign, Temp, weights,
        tab[:κ], tab[:src], scat;
        corr_χ=corr_χ, corr_S=corr_S
    )

    @. κ_ross = ifelse(χRBox > 1e-30, 1.0 / χRBox, 0.0)

    logg_val = if !isnothing(logg)
        logg
    elseif !isnothing(transition_model)
        hasproperty(transition_model, :logg) ? transition_model.logg : transition_model
    else
        nothing
    end

    wthin, wthick = if isnothing(logg_val)
        @warn "No logg was provided to the binning. The binning will only be valid for the sun!"
        wthin_l = exp.(T_type(-1.5e7) .* κ_ross)
        wthin_l, T_type(1.0) .- wthin_l
    else
        wthin_l = zeros(T_type, size(κ_ross)...)
        _thin_weight!(wthin_l, κ_ross, tab[:lnPg], ρ, T_type(exp10(logg_val)))
        wthin_l, T_type(1.0) .- wthin_l
    end

    @. opacity_table = wthin .* κBox .+ wthick .* κ_ross
    @. ϵ_table = ifelse(χBox > 1e-30, κBox / χBox, 1.0)
    @. S_table = max(SBox, T_type(1e-30))
    @. opacity_table = max(opacity_table, T_type(1e-30))

    binned = OpacityData(
        Dict{Symbol, AbstractArray}(:κ => opacity_table, :src => S_table, :ϵ => ϵ_table, :wthin => wthin),
        collect(T_type, 1:radBins); binned=true, optical_depth=false
    )

    EquationOfState(tab; opacity=binned, scattering=nothing,
                    meta=merge(filter(kv -> kv.first != :path, tab.meta), Dict{Symbol, Any}(:kind => "binned")))
end

function _thin_weight!(wthin, κ_ross, lnPg, ρ, g)
    @inbounds for k in axes(wthin, 3)
        for j in axes(wthin, 2)
            for i in axes(wthin, 1)
                τ = κ_ross[i, j, k] * exp(lnPg[i, j]) / (ρ[j] * g)
                wthin[i, j, k] = exp(eltype(wthin)(-2.0) * τ)
            end
        end
    end
    wthin
end

binning(assignment::AbstractVector{<:Integer}, weights, tab::EquationOfState; kwargs...) =
    binning(to_weight_matrix(assignment), weights, tab; kwargs...)

# ============================================================================
# Binning for 1D atmospheres
# ============================================================================

# ============================================================================
# Pre-planned 1D binning
# ============================================================================

struct Binning1DPlan{T<:AbstractFloat}
    Aχ::Matrix{T}
    AχR::Matrix{T}
    AS::Matrix{T}
    AB::Matrix{T}
    AdB::Matrix{T}
    Aκ::Matrix{T}
    ρ::Vector{T}
    Pg::Vector{T}
end

Base.size(p::Binning1DPlan) = size(p.Aχ)

function binning_1d_plan(weights, λ, ρ_1d, Temp_1d, Pg_1d, κ_1d, src_1d;
                         κ_scat_1d=nothing, corr_χ_1d=nothing, corr_S_1d=nothing)
    T_type  = eltype(κ_1d)
    Nz      = length(ρ_1d)
    nlambda = length(λ)

    Aχ  = Matrix{T_type}(undef, Nz, nlambda)
    AχR = Matrix{T_type}(undef, Nz, nlambda)
    AS  = Matrix{T_type}(undef, Nz, nlambda)
    AB  = Matrix{T_type}(undef, Nz, nlambda)
    AdB = Matrix{T_type}(undef, Nz, nlambda)
    Aκ  = Matrix{T_type}(undef, Nz, nlambda)

    @inbounds for k in 1:nlambda
        w_k = T_type(weights[k])
        λ_k = λ[k]

        @fastmath @simd for j in 1:Nz
            bnC  = T_type(Bλ_fast(λ_k, Temp_1d[j]))
            dbnC = T_type(dBdTλ_fast(λ_k, Temp_1d[j]))

            c_χ = get_corr_val(corr_χ_1d, T_type, j, k)
            c_S = get_corr_val(corr_S_1d, T_type, j, k)

            κ_val = T_type(κ_1d[j, k]) * c_χ
            inv_κ = ifelse(κ_val > 1e-30, T_type(1.0) / κ_val, T_type(0.0))
            k_abs = T_type(get_abs_k(κ_1d, κ_scat_1d, j, k)) * c_χ

            w_bnC  = w_k * bnC
            w_dbnC = w_k * dbnC

            Aχ[j, k]  = κ_val * w_bnC
            AχR[j, k] = inv_κ * w_dbnC
            AS[j, k]  = w_k * T_type(src_1d[j, k]) * c_S
            AB[j, k]  = w_bnC
            AdB[j, k] = w_dbnC
            Aκ[j, k]  = k_abs * w_bnC
        end
    end

    Binning1DPlan(Aχ, AχR, AS, AB, AdB, Aκ,
                  Vector{T_type}(T_type.(ρ_1d)), Vector{T_type}(T_type.(Pg_1d)))
end

# ============================================================================
# Submission and collection of binning tasks
# ============================================================================

function _binning_X!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
                              ρ, λ, W_assign::AbstractMatrix, Temp, weights, κ, src, κ_scat; corr_χ=nothing, corr_S=nothing)
    rhoBins = size(χBox, 2)
    AxBins  = size(χBox, 1)
    radBins = size(W_assign, 2)

    chunks = Iterators.partition(eachindex(λ), USE_BINNING_THREADS[] ? length(λ) ÷ Threads.nthreads() : length(λ)) |> collect

    begin
        tasks = map(chunks) do chunk
            if USE_BINNING_THREADS[]
                Threads.@spawn _run_bin_core!(chunk, rhoBins, AxBins, radBins, W_assign, weights, κ, κ_scat, src, Temp, λ; corr_χ=corr_χ, corr_S=corr_S)
            else
                _run_bin_core!(chunk, rhoBins, AxBins, radBins, W_assign, weights, κ, κ_scat, src, Temp, λ; corr_χ=corr_χ, corr_S=corr_S)
            end
        end
        results = fetch.(tasks)

        χBoxChunk  = getindex.(results, 1)
        χRBoxChunk = getindex.(results, 2)
        χthinChunk = getindex.(results, 3)
        SBoxChunk  = getindex.(results, 4)
        BChunk     = getindex.(results, 5)
        δBChunk    = getindex.(results, 6)
        κBoxChunk  = getindex.(results, 7)
    end

    χBox   .= sum(χBoxChunk)
    χRBox  .= sum(χRBoxChunk)
    χ_thin .= sum(χthinChunk)
    SBox   .= sum(SBoxChunk)
    B      .= sum(BChunk)
    δB     .= sum(δBChunk)
    κBox   .= sum(κBoxChunk)

    @. χBox   = ifelse(B > 1e-30, χBox / B, 0.0)
    @. κBox   = ifelse(B > 1e-30, κBox / B, 0.0)
    @. χRBox  = ifelse(δB > 1e-30, χRBox / δB, 0.0)
    @. χ_thin = ifelse(δB > 1e-30, χ_thin / δB, 0.0)

    @inbounds for j in 1:rhoBins
        χBox[:,  j, :] .*= ρ[j]
        κBox[:,  j, :] .*= ρ[j]
        χRBox[:, j, :] ./= ρ[j]
    end
end

function _run_bin_core!(chunk, rhoBins, AxBins, radBins, W_assign::AbstractMatrix, weights, κ, κ_scat, src, Temp, λ; corr_χ=nothing, corr_S=nothing)
    dbnC = Ref(eltype(κ)(0.0))
    bnC  = Ref(eltype(κ)(0.0))

    χBoxChunk  = zeros(eltype(κ), AxBins, rhoBins, radBins)
    χRBoxChunk = deepcopy(χBoxChunk)
    χthinChunk = deepcopy(χBoxChunk)
    SBoxChunk  = deepcopy(χBoxChunk)
    BChunk     = deepcopy(χBoxChunk)
    δBChunk    = deepcopy(χBoxChunk)
    κBoxChunk  = deepcopy(χBoxChunk)

    _bin_core!(
        W_assign, chunk, rhoBins, AxBins, weights, κ, κ_scat, src, Temp, λ,
        bnC, dbnC,
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

# ============================================================================
# Core binning
# ============================================================================

@inline get_corr_val(x::Nothing, T_type, I::Vararg{Int}) = T_type(1.0)
@inline get_corr_val(x::AbstractArray, T_type, I::Vararg{Int}) = T_type(x[I...])
@inline get_abs_k(κ, κ_scat::Nothing, I::Vararg{Int}) = κ[I...]
@inline get_abs_k(κ, κ_scat::AbstractArray, I::Vararg{Int}) = max(κ[I...] - κ_scat[I...], zero(eltype(κ)))

function _bin_core!(W_assign::AbstractMatrix, chunk,
                    rhoBins, AxBins, weights, κ, κ_scat, src, Temp, λ, bnC, dbnC,
                    χBoxChunk, χRBoxChunk, χthinChunk, SBoxChunk, BChunk, δBChunk, κBoxChunk,
                    corr_χ, corr_S)

    radBins = size(W_assign, 2)
    T_type = eltype(κ)

    bnC_buf  = Matrix{T_type}(undef, AxBins, rhoBins)
    dbnC_buf = Matrix{T_type}(undef, AxBins, rhoBins)

    @inbounds for k in chunk
        w_k = T_type(weights[k])
        λ_k = λ[k]

        for j in 1:rhoBins
            @simd for i in 1:AxBins
                bnC_buf[i, j]  = T_type(Bλ_fast(λ_k, Temp[i, j]))
                dbnC_buf[i, j] = T_type(dBdTλ_fast(λ_k, Temp[i, j]))
            end
        end

        for b in 1:radBins
            frac = T_type(W_assign[k, b])

            w_eff = w_k * frac

            for j in 1:rhoBins
                @simd for i in 1:AxBins
                    b_val  = bnC_buf[i, j]
                    db_val = dbnC_buf[i, j]

                    c_χ = get_corr_val(corr_χ, T_type, i, j, k)
                    c_S = get_corr_val(corr_S, T_type, i, j, k)

                    κ_val = T_type(κ[i, j, k]) * c_χ
                    s_val = T_type(src[i, j, k]) * c_S

                    w_b  = w_eff * b_val
                    w_db = w_eff * db_val

                    inv_κ = ifelse(κ_val > 1e-30, 1.0 / κ_val, T_type(0.0))

                    χBoxChunk[i, j, b]  += κ_val * w_b
                    χRBoxChunk[i, j, b] += inv_κ * w_db
                    SBoxChunk[i, j, b]  += w_eff * s_val
                    BChunk[i, j, b]     += w_b
                    δBChunk[i, j, b]    += w_db

                    k_abs = T_type(get_abs_k(κ, κ_scat, i, j, k)) * c_χ
                    κBoxChunk[i, j, b]  += k_abs * w_b
                end
            end
        end
    end

    χthinChunk .= χBoxChunk
end

# ============================================================================
# Helper functions for weight matrix conversion and opacity transition
# ============================================================================

function to_weight_matrix(assignment::AbstractVector{<:Integer})
    nbins = maximum(assignment)
    W = zeros(Float64, length(assignment), nbins)
    for i in eachindex(assignment)
        if assignment[i] > 0
            W[i, assignment[i]] = 1.0
        end
    end
    return W
end

