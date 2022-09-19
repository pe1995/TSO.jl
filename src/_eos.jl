
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
    nEiBin  = size(t.lnEi, 1)
    nRhoBin = size(t.lnRho, 2)

    lnEi2  = range(lnEi2min, lnEi2min, length=nEiBin)
    dlnEi2 = lnEi2[2] - lnEi2[1]

    lnT   = t.lnT[1, :, 1]
    lnRho = t.lnRho[:, 1, 1]

    lneo        = zeros(nEiBin)
    lnT2        = zeros(nEiBin)
    dThdEi2     = zeros(nEiBin)
    dlnEdlnRho2 = zeros(nEiBin)
    
    Tg2     = zeros(nEiBin, nRhoBin)
    lnNe2   = zeros(nEiBin, nRhoBin)
    lnRoss2 = zeros(nEiBin, nRhoBin)
    lnPg2   = zeros(nEiBin, nRhoBin)

    f    = zeros(nEiBin)
    g    = zeros(nEiBin)
    d    = zeros(nEiBin)
    ff   = zeros(nEiBin)
    dd   = zeros(nEiBin)
    gg   = zeros(nEiBin)
    hh   = zeros(nEiBin)

    for i in eachindex(lnRho)
        etmin = 1.0001*minimum(view(t.lnEi, :, i, 1))
        etmax = 0.9999*maximum(view(t.lnEi, :, i, 1))

        # index range of the new axis to consider
        esi = findfirst(x->x>etmin, lnEi2)
        eei = findfirst(x->x>etmax, lnEi2)

        (isnothing(esi) | isnothing(eei)) && error("empty ei bin at given rho") 

        # first do the ln Tg -> ln(ei) swap ###################################
        lneo .= t.lnEi[:, i, 1]         # ln ei on lnT grid
        f    .= lnT                     # lnT grid values
        d    .= 1.0 ./ t.lnEi[:, i, 3]  # d lnT / d ln(ei) | rho on lnT grid
        g    .= t.lnEi[:, i, 2]         # d ln ei / d ln(rho)| Tg on lnT grid
       
        # get h = d/d ln ei ( d ln ei / d ln rho) on lnT grid
        s   = spline(interpolate(lneo, g, BSplineOrder(4)))
        h   = (Derivertive(1) * s).(lneo)
        
        # Interpolate to the new grid
        #ff .= linear_interpolation(lneo, f, extrapolation_bc=Line()).(lnEi2)
        #dd .= linear_interpolation(lneo, d, extrapolation_bc=Line()).(lnEi2)
        #gg .= linear_interpolation(lneo, g, extrapolation_bc=Line()).(lnEi2)
        #hh .= linear_interpolation(lneo, h, extrapolation_bc=Line()).(lnEi2)

        # Tabgen interpolates the above using the spline derivative, dont know why though
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnEi2, lneo, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnEi2, lneo, g, h)

        lnT2        .= ff    # lnT on ln(ei) grid
        dThdEi2     .= dd    # d lnT /d ln (ei) on ln(ei) grid
        dlnEdlnRho2 .= gg    # d ln ei / d ln(rho) | Tg   on ln(ei) grid

        ff .= exp.(lnT2)     # Tg on ln(ei) grid
        dd .= -ff            # d Tg d lnT on ln(ei) grid
        gg .= 0.0
        Tg2[:, i] .= ff

        add_outside!(view(Tg2, :, i), esi, eei)

        # Now the pressure ####################################################
        f .= t.lnPg[:, i, 1] # f on lnT grid
        d .= t.lnPg[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnPg[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivertive(1) * s).(lnT)

        #ff .= linear_interpolation(lnT, f, extrapolation_bc=Line()).(lnT2)
        #dd .= linear_interpolation(lnT, d, extrapolation_bc=Line()).(lnT2)
        #gg .= linear_interpolation(lnT, g, extrapolation_bc=Line()).(lnT2)
        #hh .= linear_interpolation(lnT, h, extrapolation_bc=Line()).(lnT2)
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

        lnPg2[:, i] .= ff 
        add_outside!(view(lnPg2, :, i), esi, eei)

        # the electron density ################################################
        f .= t.lnNe[:, i, 1] # f on lnT grid
        d .= t.lnNe[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnNe[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivertive(1) * s).(lnT)
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

        lnNe2[:, i] .= ff 
        add_outside!(view(lnNe2, :, i), esi, eei)

        # the rosseland opacity ###############################################
        f .= t.lnRoss[:, i, 1] # f on lnT grid
        d .= t.lnRoss[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnRoss[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivertive(1) * s).(lnT)
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

        lnRoss2[:, i] .= ff 
        add_outside!(view(lnRoss2, :, i), esi, eei)
    end

    RegularEoSTable(lnRho, lnT2, lnEi2, lnPg2, lnRoss2, lnNe2)
end

lnT_to_lnEi(t) = error("No transformation from T to Ei implemented for this table.")

"""Interpolate the function and first derivative of a cubic spline. (Taken from Tabgen)"""
function interpolate_f_df!(ff, dd, xx, x, f, d)
    @inbounds for i in eachindex(ff)
        if (x[2]>x[1])
            for k=2:n
                (xx[i] < x[k]) && break
            end
        else
            for k=2:n
                (x[k] < xx[i]) && break
            end
        end

        @fastmath begin 
            dx  = x[k] - x[k-1]
            p   = (xx[i] - x[k-1])/dx
            q   = 1. - p
            qp  = q*p
            qmp = q - p
            pf  = p - qp*qmp
            pd  = -p*qp
            qd  = q*qp
            pdf = 1. -qmp*qmp + 2.0*qp
            pdd = p*(-qmp - q)
            qdd = q*(qmp - p)
            ff[i] = f[k-1] + pf*(f[k]-f[k-1]) + qd*d[k-1]*dx + pd*d[k]*dx
            dd[i] = pdf*(f[k]-f[k-1])/dx + qdd*d[k-1] + pdd*d[k]
        end
    end
end

"""Add linear extrapolation outside of the given index."""
function add_outside!(f::Vector{K}, esi, eei) where {K}
    nF = length(f)
    for k=1:esi-1
        f[esi-k] = 2.0*f[esi] - f[esi+k]
    end 
    
    for k=1:nF-eei
        f[eei+k] = 2.0*f[eei] - f[eei-k]
    end 
end