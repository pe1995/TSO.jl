
"""Switch between different grids of the EoS table."""
function transform(t::EoSTable, from_to::Pair{Symbol, Symbol})
    return if (first(from_to) == :lnT) & (last(from_to) == :lnEi)
        lnT_to_lnEi(t)
    else
        error("The given transformation of grids is not available yet.")
        nothing
    end
end

"""
Switch from a temperature to an internal energy grid.
This routine will invert the table, going from ln rho and lnT as
independent variables to ln rho and ln ei as independent variables.

This is a bit tricky because d f / d ln(rho) at constant lnT is
not equal to d f / d ln(rho) at constant ln(ei).

To make everything a function of ln(E) instead of lnT, use the relation

                d ln(E)          d ln(E)
    d ln(E) =  ------- dlnT + ---------  d ln(rho) = 0   when
                dlnT           d ln(rho)

        dlnT       (d ln(E)/d ln(rho))
    --------- = - -------------------
    d ln(rho)      (d ln(E)/dlnT)

Thus

        df                            df                      df
    ------ at constant ln(E) is = ------ at const lnT - ------ * above
    dlnrho                        dlnrho                  dlnT

esi,eei contain the Ei index below and above which the
temperature is outside the range in the ln rho-lnT table for a
given rho(i).
"""
function lnT_to_lnEi(t::RegularEoSTable)
    # new axis range
    lnEi2min = minimum(t.lnEi)
    lnEi2max = maximum(t.lnEi)
    lnEi2r   = lnEi2max - lnEi2min

    # new regular axis
    nEiBin = size(t.lnEi, 1)
    lnEi2  = range(lnEi2min, lnEi2min, length=nEiBin)
    dlnEi2 = lnEi2[2] - lnEi2[1]

    lnT   = t.lnT[1, :, 1]
    lnRho = t.lnRho[:, 1, 1]

    lneo = zeros(nEiBin)
    f    = zeros(nEiBin)
    g    = zeros(nEiBin)
    d    = zeros(nEiBin)

    for i in eachindex(t.lnRho)
        etmin = 1.0001*minimum(view(t.lnEi, :, i, 1))
        etmax = 0.9999*maximum(view(t.lnEi, :, i, 1))

        # index range of the new axis to consider
        esi = findfirst(x->x>etmin, lnEi2)
        eei = findfirst(x->x>etmax, lnEi2)

        (isnothing(esi) | isnothing(eei)) && error("empty ei bin at given rho") 

        # first do the ln Tg -> ln(ei) swap
        lneo .= t.lnEi[:, i, 1]         # ln ei on lnT grid
        f    .= lnT                     # lnT grid values
        d    .= 1.0 ./ t.lnEi[:, i, 3]  # d lnT / d ln(ei) | rho on lnT grid
        g    .= t.lnEi[:, i, 2]         # d ln ei / d ln(rho)| Tg on lnT grid
       
        # get h = d/d ln ei ( d ln ei / d ln rho) on lnT grid
        s   = spline(interpolate(lneo, g, BSplineOrder(4)))
        h   = (Derivertive(1) * s).(lneo)
        

    end
end

lnT_to_lnEi(t) = error("No transformation from T to Ei implemented for this table.")