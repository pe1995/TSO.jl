# ============================================================================
# Clean data-structure variants for T -> Ei conversion
# ============================================================================

"""
    remap_T_to_E(tab; upsample=1000, extend=false, lnEimin=nothing, lnEimax=nothing)

Remap a table from the temperature grid to the internal energy grid, as needed by DISPATCH.
"""
function remap_T_to_E(tab::EquationOfState;
    upsample=1000, extend=false,
    downD=0.4, downE=0.1, upD=0.01, upE=0.01,
    lnEimin=nothing, lnEimax=nothing)

    tab_new = if extend
        TSO.extend(tab, downD=downD, downE=downE, upD=upD, upE=upE)
    else
        EquationOfState(tab; eos=EoSData(Dict{Symbol, AbstractArray}(k => copy(v) for (k, v) in tab.eos.data)))
    end

    smoothAccumulate!(tab_new, spline=true)

    if !isnothing(lnEimax)
        lnEi = tab_new[:lnEi]
        lnEi[lnEi .>= lnEimax] .= lnEimax
        isnothing(lnEimin) || (lnEi[lnEi .<= lnEimin] .= lnEimin)
    end

    tabE = switch_energy(tab_new, upsample=upsample)

    for p in (:opacity, :scattering)
        oE, o = getfield(tabE, p), getfield(tab_new, p)
        isnothing(oE) && continue
        for f in _float_fields(oE), i in axes(oE[f], 3)
            ndims(oE[f]) == 3 || continue
            limit_beyond!(view(oE[f], :, :, i), view(o[f], :, :, i))
        end
    end
    for f in _eos_fields(tabE)
        haskey(tab_new, f) && limit_beyond!(tabE.eos[f], tab_new[f])
    end

    fill_nan!(tabE)

    return tabE
end

function get_e_limit(eos::EquationOfState, model_lnEi::AbstractArray, lnEi_stretch=1.0)
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
