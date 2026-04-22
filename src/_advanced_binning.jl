# ============================================================================
# Binning
# ============================================================================

const USE_BINNING_THREADS = Ref(true)

function advanced_binning(W_assign::AbstractMatrix, weights, aos::E, opacities, scattering=nothing; 
                          transition_model=nothing, logg=nothing, corr_χ=nothing, corr_S=nothing) where {E<:AxedEoS}
    
    radBins = size(W_assign, 2)
    eos = aos.eos
    eaxis = is_internal_energy(aos)

    iscat = isnothing(scattering)
    remove_from_thin = !iscat

    if !remove_from_thin
        @warn "No scattering opacity provided!"
    end

    rhoBins = length(eos.lnRho)
    AxBins  = aos.energy_axes.length
    T_type  = eltype(eos.lnRho)

    # Allocate required arrays locally
    χBox   = zeros(T_type, AxBins, rhoBins, radBins)
    χRBox  = zeros(T_type, AxBins, rhoBins, radBins)
    κ_ross = zeros(T_type, AxBins, rhoBins, radBins)
    χ_thin = zeros(T_type, AxBins, rhoBins, radBins)
    κBox   = zeros(T_type, AxBins, rhoBins, radBins)
    SBox   = zeros(T_type, AxBins, rhoBins, radBins)
    B      = zeros(T_type, AxBins, rhoBins, radBins)
    δB     = zeros(T_type, AxBins, rhoBins, radBins)

    Temp = zeros(T_type, AxBins, rhoBins)
    if ndims(eos.lnT) == 2
        Temp .= exp.(eos.lnT)
    else
        puff_up!(Temp, exp.(eos.lnT))
    end

    ρ = exp.(eos.lnRho)

    #if Threads.nthreads() > 1
    #    @info "Binning opacities with $(Threads.nthreads()) threads."
    #end

    scat = remove_from_thin ? scattering.κ : nothing
    
    # Run the multithreaded orchestrator
    _advanced_binning_X!(
        B, δB, SBox, κBox, χBox, χRBox, χ_thin, ρ,
        opacities.λ, W_assign, Temp, weights, 
        opacities.κ, opacities.src, scat;
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
        @warn "No logg was provided to the binning. " *
              "The binning will only be valid for the sun!"
        wthin_l = exp.(T_type(-1.5e7) .* κ_ross)
        wthin_l, T_type(1.0) .- wthin_l
    else
        g = T_type(exp10(logg_val))
        wthin_l = zeros(T_type, size(κ_ross)...)
        @inbounds for k in 1:radBins
            for j in 1:rhoBins
                for i in 1:AxBins
                    τ = κ_ross[i, j, k] * exp(eos.lnPg[i, j]) / (ρ[j] * g)
                    wthin_l[i, j, k] = exp(T_type(-2.0) * τ)
                end
            end
        end
        wthin_l, T_type(1.0) .- wthin_l
    end

    opacity_table = wthin .* κBox .+ wthick .* κ_ross     
    ϵ_table = @. ifelse(χBox > 1e-30, κBox / χBox, 1.0)
    S_table = SBox 

    return BinnedOpacities(
        SqOpacity(
            opacity_table, opacities.κ_ross, S_table,
            collect(T_type, 1:radBins), false
        ), 
        ϵ_table, 
        wthin
    )
end

function advanced_binning(binning::AbstractVector{<:Integer}, weights, aos::E, opacities, scattering=nothing; kwargs...) where {E<:AxedEoS}
    W_assign = to_weight_matrix(binning)
    return advanced_binning(W_assign, weights, aos, opacities, scattering; kwargs...)
end

# ============================================================================
# Binning for 1D atmospheres
# ============================================================================

function advanced_binning_1d(W_assign::AbstractMatrix, weights, λ, ρ_1d, Temp_1d, Pg_1d, κ_1d, src_1d; 
                             κ_scat_1d=nothing, corr_χ_1d=nothing, corr_S_1d=nothing, logg=nothing)
    
    radBins = size(W_assign, 2)
    Nz = length(ρ_1d)
    nlambda = length(λ)
    T_type = eltype(κ_1d)
    
    # Allocate locally
    χBox   = zeros(T_type, 1, Nz, radBins)
    χRBox  = zeros(T_type, 1, Nz, radBins)
    κ_ross = zeros(T_type, 1, Nz, radBins)
    χ_thin = zeros(T_type, 1, Nz, radBins)
    κBox   = zeros(T_type, 1, Nz, radBins)
    SBox   = zeros(T_type, 1, Nz, radBins)
    B      = zeros(T_type, 1, Nz, radBins)
    δB     = zeros(T_type, 1, Nz, radBins)

    Temp_2d = reshape(Temp_1d, 1, Nz)
    κ_3d    = reshape(κ_1d, 1, Nz, nlambda)
    src_3d  = reshape(src_1d, 1, Nz, nlambda)
    scat_3d = isnothing(κ_scat_1d) ? nothing : reshape(κ_scat_1d, 1, Nz, nlambda)
    corr_X  = isnothing(corr_χ_1d) ? nothing : reshape(corr_χ_1d, 1, Nz, nlambda)
    corr_S  = isnothing(corr_S_1d) ? nothing : reshape(corr_S_1d, 1, Nz, nlambda)
    
    _advanced_binning_X!(
        B, δB, SBox, κBox, χBox, χRBox, χ_thin,
        ρ_1d, λ, W_assign, Temp_2d, weights, κ_3d, src_3d, scat_3d; 
        corr_χ=corr_X, corr_S=corr_S
    )
    
    @. κ_ross = ifelse(χRBox > 1e-30, 1.0 / χRBox, 0.0)
    
    wthin, wthick = if isnothing(logg)
        wthin_l = exp.(T_type(-1.5e7) .* κ_ross)
        wthin_l, T_type(1.0) .- wthin_l
    else
        opacity_transition_1d(κ_ross, ρ_1d, Pg_1d, logg)
    end
    
    @. κBox = wthin * κBox + wthick * κ_ross 
    
    # Returns the updated binned arrays directly
    return (κBox=κBox, SBox=SBox, χBox=χBox, χRBox=χRBox)
end

function advanced_binning_1d_quick(W_assign::AbstractMatrix, weights, λ, ρ_1d, Temp_1d, Pg_1d, κ_1d, src_1d; 
                             κ_scat_1d=nothing, corr_χ_1d=nothing, corr_S_1d=nothing, logg=nothing)
    
    radBins = size(W_assign, 2)
    Nz      = length(ρ_1d)
    nlambda = length(λ)
    T_type  = eltype(κ_1d)
    
    χBox  = zeros(T_type, 1, Nz, radBins)
    χRBox = zeros(T_type, 1, Nz, radBins)
    κBox  = zeros(T_type, 1, Nz, radBins)
    SBox  = zeros(T_type, 1, Nz, radBins)
    
    B_norm  = zeros(T_type, 1, Nz, radBins)
    dB_norm = zeros(T_type, 1, Nz, radBins)

    bnC_col   = Vector{T_type}(undef, Nz)
    dbnC_col  = Vector{T_type}(undef, Nz)
    κ_col     = Vector{T_type}(undef, Nz)
    inv_κ_col = Vector{T_type}(undef, Nz)
    src_col   = Vector{T_type}(undef, Nz)
    κ_abs_col = Vector{T_type}(undef, Nz)

    @inbounds for k in 1:nlambda
        w_k = T_type(weights[k])
        λ_k = λ[k]
        
        @inbounds @fastmath @simd for j in 1:Nz
            bnC_col[j]  = T_type(Bλ_fast(λ_k, Temp_1d[j]))
            dbnC_col[j] = T_type(dBdTλ_fast(λ_k, Temp_1d[j]))
            
            c_χ = get_corr_val(corr_χ_1d, T_type, j, k)
            c_S = get_corr_val(corr_S_1d, T_type, j, k)

            κ_val = T_type(κ_1d[j, k]) * c_χ
            κ_col[j] = κ_val
            inv_κ_col[j] = ifelse(κ_val > 1e-30, 1.0 / κ_val, T_type(0.0))
            src_col[j] = T_type(src_1d[j, k]) * c_S
            
            k_abs = T_type(get_abs_k(κ_1d, κ_scat_1d, j, k))
            κ_abs_col[j] = k_abs * c_χ
        end
        
        for b in 1:radBins
            frac = T_type(W_assign[k, b])
            #(frac <= 1e-8) && continue 
            
            w_eff = w_k * frac
            
            @turbo for j in 1:Nz
                w_bnC   = w_eff * bnC_col[j]
                w_dbnC  = w_eff * dbnC_col[j]
                
                χBox[1, j, b]    += κ_col[j] * w_bnC
                χRBox[1, j, b]   += inv_κ_col[j] * w_dbnC
                SBox[1, j, b]    += w_eff * src_col[j]
                
                B_norm[1, j, b]  += w_bnC
                dB_norm[1, j, b] += w_dbnC
                
                κBox[1, j, b]    += κ_abs_col[j] * w_bnC
            end
        end
    end

    @inbounds for b in 1:radBins
        @turbo for j in 1:Nz
            ρ_j = T_type(ρ_1d[j])
            bn  = B_norm[1, j, b]
            dbn = dB_norm[1, j, b]

            χBox[1, j, b] = ifelse(bn > 1e-30, (χBox[1, j, b] / bn) * ρ_j, T_type(0.0))
            κBox[1, j, b] = ifelse(bn > 1e-30, (κBox[1, j, b] / bn) * ρ_j, T_type(0.0))
            
            χR_tmp = ifelse(dbn > 1e-30, χRBox[1, j, b] / dbn, T_type(0.0))
            χRBox[1, j, b] = ifelse(χR_tmp > 1e-30, χR_tmp / ρ_j, T_type(0.0))
        end
    end

    if isnothing(logg)
        @inbounds @simd for i in eachindex(κBox)
            χR_val = χRBox[i]
            κ_ross = ifelse(χR_val > 1e-30, 1.0 / χR_val, T_type(0.0))
            
            wthin  = exp(T_type(-1.5e7) * κ_ross)
            κBox[i] = wthin * κBox[i] + (T_type(1.0) - wthin) * κ_ross
        end
    else
        g = T_type(exp10(logg))
        
        @inbounds for b in 1:radBins
            @simd for z in 1:Nz
                χR_val = χRBox[1, z, b]
                κ_ross = ifelse(χR_val > 1e-30, 1.0 / χR_val, T_type(0.0))
                
                Pg_val = T_type(Pg_1d[z])
                ρ_val  = T_type(ρ_1d[z])
                τ_val  = (κ_ross * Pg_val) / (ρ_val * g)
                
                wthin  = exp(T_type(-2.0) * τ_val)
                wthick = T_type(1.0) - wthin
                
                κBox[1, z, b] = wthin * κBox[1, z, b] + wthick * κ_ross
            end
        end
    end

    return (κBox=κBox, SBox=SBox, χBox=χBox, χRBox=χRBox)
end

# ============================================================================
# Submission and collection of binning tasks
# ============================================================================

function _advanced_binning_X!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
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

    # Sum over all chunks into the top-level arrays
    χBox   .= sum(χBoxChunk)
    χRBox  .= sum(χRBoxChunk)
    χ_thin .= sum(χthinChunk) 
    SBox   .= sum(SBoxChunk)
    B      .= sum(BChunk)
    δB     .= sum(δBChunk)
    κBox   .= sum(κBoxChunk)

    # Prevent 0.0 / 0.0 = NaN for empty bins using safe division
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
@inline get_abs_k(κ, κ_scat::AbstractArray, I::Vararg{Int}) = κ[I...] - κ_scat[I...]

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
            #(frac <= 1e-8) && continue 
            
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

                    inv_κ = 1.0 / κ_val
                    
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

function to_weight_matrix(binning::AbstractVector{<:Integer})
    nbins = maximum(binning)
    W = zeros(Float64, length(binning), nbins)
    for i in eachindex(binning)
        if binning[i] > 0
            W[i, binning[i]] = 1.0
        end
    end
    return W
end

function opacity_transition_1d(κ_ross_1d, ρ_1d, Pg_1d, logg)
    T_type = eltype(κ_ross_1d)
    
    # κ_ross_1d is typically a view or 3D array sliced.
    _, Nz, radBins = size(κ_ross_1d)
    
    τ = zeros(T_type, 1, Nz, radBins)
    g = exp10(logg)
    
    @inbounds for b in 1:radBins
        for z in 1:Nz
            # τ = κ * P_g / (ρ * g)
            τ[1, z, b] = κ_ross_1d[1, z, b] * Pg_1d[z] / (ρ_1d[z] * g)
        end
    end

    wthin = exp.(T_type(-2.0) .* τ)
    return wthin, T_type(1.0) .- wthin
end