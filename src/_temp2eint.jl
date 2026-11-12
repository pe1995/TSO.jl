# ============================================================================
# Clean data-structure variants for T -> Ei conversion
# ============================================================================

function remap_T_to_E(aos::AxedEoS, opa::SqOpacity;
    upsample=1000, extend=false,
    downD=0.4, downE=0.1, upD=0.01, upE=0.01,
    lnEimin=nothing, lnEimax=nothing)

    eos_new, opa_new = if extend
        TSO.extend(aos, opa, downD=downD, downE=downE, upD=upD, upE=upE)
    else
        aos.eos, opa
    end

    TSO.smoothAccumulate!(@axed(eos_new), spline=true)

    if !isnothing(lnEimax)
        lnEimax_orig = maximum(eos_new.lnEi)
        mask = eos_new.lnEi.>=lnEimax
        eos_new.lnEi[mask] .= lnEimax

        if !isnothing(lnEimin)
            eos_new.lnEi[eos_new.lnEi.<=lnEimin] .= lnEimin
        end
    end

    eosE, opaE = switch_energy(@axed(eos_new), opa_new, upsample=upsample, conservative=false)
    aosE = @axed eosE

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

function get_e_limit(eos::SqEoS, model_lnEi::AbstractArray, lnEi_stretch=1.0)
    lnEimin = minimum(model_lnEi)
    lnEimax = maximum(model_lnEi)

    tmin = minimum(eos.lnT)

    lnEi_min_curve = eos.lnEi[1, :]

    lnEimin = median(lnEi_min_curve)

    return (lnEimin, lnEimax + abs(lnEimax - lnEimin) * lnEi_stretch)
end

function limit_beyond!(new_arr, old_arr; set_to=NaN)
    m = (new_arr .> maximum(old_arr)) .|| (new_arr .< minimum(old_arr))
    new_arr[m] .= set_to
end
