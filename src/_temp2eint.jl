# ============================================================================
# Clean data-structure variants for T -> Ei conversion
# ============================================================================

"""
    remap_T_to_E(aos::AxedEoS, opa::SqOpacity; kwargs...)

Convert the given AxedEoS and SqOpacity from a Temperature grid to an Internal Energy grid structurally, without any file I/O.
Returns the newly interpolated `(eosE, opaE)`.
"""
function remap_T_to_E(aos::AxedEoS, opa::SqOpacity;
    upsample=1000, extend=false, 
    downD=0.4, downE=0.1, upD=0.01, upE=0.01, 
    lnEimin=nothing, lnEimax=nothing)

    # possibly extrapolate EoS
    eos_new, opa_new = if extend
        #@info "Extrapolating EoS beyond limits."
        TSO.extend(aos, opa, downD=downD, downE=downE, upD=upD, upE=upE)
    else
        aos.eos, opa
    end

    # Again making it monotonic
    TSO.smoothAccumulate!(@axed(eos_new), spline=true)

    # Limit the internal energy range to what is occupied by the initial model to make sure
    # the model atmosphere will be well sampled
    if !isnothing(lnEimax)
        lnEimax_orig = maximum(eos_new.lnEi)
        mask = eos_new.lnEi.>=lnEimax
        eos_new.lnEi[mask] .= lnEimax

        if !isnothing(lnEimin)
            #@info "Limiting the internal energy (log) between $(lnEimin) - $(lnEimax) on the energy grid." 
            #@info "The original limits were $(minimum(eos_new.lnEi)) - $(lnEimax_orig)." 
            eos_new.lnEi[eos_new.lnEi.<=lnEimin] .= lnEimin
        else
            #@info "Limiting the internal energy (log) to $(lnEimax) on the energy grid." 
            #@info "The original limits was $(lnEimax_orig)." 
        end
        #@info "The upper limit will have an effect on $(count(mask)/length(mask)*100)% of points in the rho-T table."
    end

    eosE, opaE = switch_energy(@axed(eos_new), opa_new, upsample=upsample, conservative=false)
    aosE = @axed eosE

    # do not allow values larger than in the original table
    for (i, l) in enumerate(opaE.λ)
        limit_beyond!(view(opaE.κ, :, :, i), view(opa_new.κ, :, :, i))
        limit_beyond!(view(opaE.src, :, :, i), view(opa_new.src, :, :, i))
    end
    limit_beyond!(eosE.lnPg, eos_new.lnPg)
    limit_beyond!(eosE.lnNe, eos_new.lnNe)
    limit_beyond!(eosE.lnRoss, eos_new.lnRoss)
    limit_beyond!(eosE.lnT, eos_new.lnT)
    
    TSO.fill_nan!(aosE, opaE)
    
    return eosE, opaE
end

"""
    get_e_limit(eos::SqEoS, model_lnEi::AbstractArray, lnEi_stretch=1.0)

Extract the internal energy boundaries given an EoS and a stellar model's internal energy generically.
"""
function get_e_limit(eos::SqEoS, model_lnEi::AbstractArray, lnEi_stretch=1.0)
    lnEimin = minimum(model_lnEi)
    lnEimax = maximum(model_lnEi)

    # determin the temperature minimum
    tmin = minimum(eos.lnT)

    # lookup the corresponding internal energy curve
    lnEi_min_curve = eos.lnEi[1, :]

    # use the median to avoid overly small lower limits
    lnEimin = median(lnEi_min_curve)

    return (lnEimin, lnEimax + abs(lnEimax - lnEimin) * lnEi_stretch)
end
