#========================== functionality ====================================#

"""
Switch between different grids of the EoS table.
"""
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
function lnT_to_lnEi(t::RegularEoSTable, interpolate_with_derivative=true)
    T = eltype(t.lnEi)
    # new axis range
    lnEi2min = minimum(view(t.lnEi, :, :, 1))
    lnEi2max = maximum(view(t.lnEi, :, :, 1))
    lnEi2r   = lnEi2max - lnEi2min

    # new regular axis
    nEiBin  = size(t.lnEi, 1)
    nRhoBin = length(t.lnRho)

    lnEi2  = Base.convert.(T, range(lnEi2min, lnEi2max, length=nEiBin))
    dlnEi2 = lnEi2[2] - lnEi2[1]

    lnT   = t.lnT
    lnRho = t.lnRho

    mask        = zeros(Int, nEiBin)
    lneo        = zeros(T, nEiBin)
    lnT2        = zeros(T, nEiBin)
    dThdEi2     = zeros(T, nEiBin)
    dlnEdlnRho2 = zeros(T, nEiBin)
    
    Tg2     = zeros(T, nEiBin, nRhoBin)
    lnNe2   = zeros(T, nEiBin, nRhoBin)
    lnRoss2 = zeros(T, nEiBin, nRhoBin)
    lnPg2   = zeros(T, nEiBin, nRhoBin)

    f    = zeros(T, nEiBin)
    g    = zeros(T, nEiBin)
    d    = zeros(T, nEiBin)
    ff   = zeros(T, nEiBin)
    dd   = zeros(T, nEiBin)
    gg   = zeros(T, nEiBin)
    hh   = zeros(T, nEiBin)

    esi ::Union{Int,Nothing} = nothing
    eei ::Union{Int,Nothing} = nothing

    for i in eachindex(lnRho)
        etmin = 1.0001*minimum(view(t.lnEi, :, i, 1))
        etmax = 0.9999*maximum(view(t.lnEi, :, i, 1))

        # index range of the new axis to consider
        esi = findfirst(x->x>etmin, lnEi2)
        eei = findfirst(x->x>etmax, lnEi2)

        (isnothing(esi) | isnothing(eei)) && error("empty ei bin at given rho") 

        # first do the ln Tg -> ln(ei) swap ###################################
        mask .= sortperm(t.lnEi[:, i, 1])
        lneo .= t.lnEi[mask, i, 1]         # ln ei on lnT grid
        f    .= lnT[mask]                  # lnT grid values
        d    .= 1.0 ./ t.lnEi[mask, i, 3]  # d lnT / d ln(ei) | rho on lnT grid
        g    .= t.lnEi[mask, i, 2]         # d ln ei / d ln(rho)| Tg on lnT grid
       
        # get h = d/d ln ei ( d ln ei / d ln rho) on lnT grid
        s   = spline(interpolate(lneo, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lneo)
        
        # Interpolate to the new grid
        #ff .= linear_interpolation(lneo, f, extrapolation_bc=Line()).(lnEi2)
        #dd .= linear_interpolation(lneo, d, extrapolation_bc=Line()).(lnEi2)
        #gg .= linear_interpolation(lneo, g, extrapolation_bc=Line()).(lnEi2)
        #hh .= linear_interpolation(lneo, h, extrapolation_bc=Line()).(lnEi2)

        # Tabgen interpolates the above using the spline derivative, dont know why though
        interpolate_f_df!(ff, dd, lnEi2, lneo, f, d)
        interpolate_f_df!(gg, hh, lnEi2, lneo, g, h)

        lnT2        .= ff    # lnT on ln(ei) grid
        dThdEi2     .= dd    # d lnT /d ln (ei) on ln(ei) grid
        dlnEdlnRho2 .= gg    # d ln ei / d ln(rho) | Tg   on ln(ei) grid

        ff .= exp.(lnT2)     # Tg on ln(ei) grid
        dd .= -ff            # d Tg d lnT on ln(ei) grid
        gg .= 0.0
        Tg2[:, i] .= ff

        if !interpolate_with_derivative
            mask = sortperm(view(lnEi2, esi:eei))
            x  = view(lnEi2, esi:eei)[mask]
            y  = view(ff, esi:eei)[mask]
            ip = linear_interpolation(x, y, extrapolation_bc=Line())
            Tg2[1:esi-1, i]   .= ip(view(lnEi2, 1:esi-1))
            Tg2[eei+1:end, i] .= ip(view(lnEi2, eei+1:length(lnEi2)))
        else
            add_outside!(view(Tg2, :, i), esi, eei)
        end

        # Now the pressure ####################################################
        f .= t.lnPg[:, i, 1] # f on lnT grid
        d .= t.lnPg[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnPg[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lnT)

        #ff .= linear_interpolation(lnT, f, extrapolation_bc=Line()).(lnT2)
        #dd .= linear_interpolation(lnT, d, extrapolation_bc=Line()).(lnT2)
        #gg .= linear_interpolation(lnT, g, extrapolation_bc=Line()).(lnT2)
        #hh .= linear_interpolation(lnT, h, extrapolation_bc=Line()).(lnT2)
        
        interpolate_f_df!(ff, dd, lnT2, lnT, f, d)
        interpolate_f_df!(gg, hh, lnT2, lnT, g, h)

        lnPg2[:, i] .= ff 

        if !interpolate_with_derivative
            mask = sortperm(view(lnT2, esi:eei))
            x = view(lnT2, esi:eei)[mask]
            y = view(ff, esi:eei)[mask]
            ip = linear_interpolation(x, y, extrapolation_bc=Line())
            lnPg2[1:esi-1, i]   .= ip.(view(lnT2, 1:esi-1))
            lnPg2[eei+1:end, i] .= ip.(view(lnT2, eei+1:length(lnT2)))
        else
            add_outside!(view(lnPg2, :, i), esi, eei)
        end

        # the electron density ################################################
        f .= t.lnNe[:, i, 1] # f on lnT grid
        d .= t.lnNe[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnNe[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lnT)
 
        interpolate_f_df!(ff, dd, lnT2, lnT, f, d)
        interpolate_f_df!(gg, hh, lnT2, lnT, g, h)

        lnNe2[:, i] .= ff 

        if !interpolate_with_derivative
            mask = sortperm(view(lnT2, esi:eei))
            x = view(lnT2, esi:eei)[mask]
            y = view(ff, esi:eei)[mask]
            ip = linear_interpolation(x, y, extrapolation_bc=Line())
            lnNe2[1:esi-1, i]   .= ip.(view(lnT2, 1:esi-1))
            lnNe2[eei+1:end, i] .= ip.(view(lnT2, eei+1:length(lnT2)))
        else
            add_outside!(view(lnNe2, :, i), esi, eei)
        end

        if any(isnan.(view(lnNe2, :, i)))
            @show i  "after" view(lnNe2, :, i)
        end

        # the rosseland opacity ###############################################
        f .= t.lnRoss[:, i, 1] # f on lnT grid
        d .= t.lnRoss[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnRoss[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lnT)
        
        interpolate_f_df!(ff, dd, lnT2, lnT, f, d)
        interpolate_f_df!(gg, hh, lnT2, lnT, g, h)

        lnRoss2[:, i] .= ff 

        if !interpolate_with_derivative
            mask = sortperm(view(lnT2, esi:eei))
            x = view(lnT2, esi:eei)[mask]
            y = view(ff, esi:eei)[mask]
            ip = linear_interpolation(x, y, extrapolation_bc=Line())
            lnRoss2[1:esi-1, i]   .= ip.(view(lnT2, 1:esi-1))
            lnRoss2[eei+1:end, i] .= ip.(view(lnT2, eei+1:length(lnT2)))
        else
            add_outside!(view(lnRoss2, :, i), esi, eei)
        end
    end

    RegularEoSTable(lnRho, log.(Tg2), lnEi2, lnPg2, lnRoss2, lnNe2)
end

nanmin(x) = isnan(x) ? Inf  : x
nanmax(x) = isnan(x) ? -Inf : x

"""
A simple method to interpolate the EoS to new Energy grid.
"""
function energy_grid(t::RegularEoSTable)
    T = eltype(t.lnEi)

    # new axis range
    #@show size(t.lnEi)
    lnEi2min = minimum(nanmin, t.lnEi)
    lnEi2max = maximum(nanmax, t.lnEi)

    #@show maximum(nanabs, minimum(t.lnEi, dims=2)) minimum(nanabs, maximum(t.lnEi, dims=2))
    #@show minimum(nanmin, t.lnEi, dims=2) maximum(nanmax, t.lnEi, dims=2)

    # The upper limit should not correspond to a temperature larger than 500,000K (total upper limit)
    #e_5 = first(lookup(t, :lnEi, [t.lnRho[end]], [log(500000)]))
    #lnEi2max = min(lnEi2max, e_5)

    # new regular axis
    nEiBin  = size(t.lnEi, 1)
    nRhoBin = length(t.lnRho)

    lnEi2  = Base.convert.(T, range(lnEi2min, lnEi2max, length=nEiBin))
    dlnEi2 = lnEi2[2] - lnEi2[1]

    # New arrays
    lnTg2   = zeros(T, nEiBin, nRhoBin)
    lnNe2   = zeros(T, nEiBin, nRhoBin)
    lnRoss2 = zeros(T, nEiBin, nRhoBin)
    lnPg2   = zeros(T, nEiBin, nRhoBin)

    # old T axis as range for cubic interpolation
    old_t = range(minimum(t.lnT), maximum(t.lnT), length=length(t.lnT))

    x = zeros(nEiBin)
    y = zeros(nEiBin)
    t_axis = zeros(nEiBin)
    mask = zeros(Int, nEiBin)

    for i in eachindex(t.lnRho)
        # This is a column of the table, interpolate it and evaluate at the right position

        # It is listed as a function of temperature -> get the new T axis for the new E axis
        mask .= sortperm(t.lnEi[:, i])
        x    .= t.lnEi[mask, i]
        y    .= t.lnT[mask]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        t_axis      .= ip.(lnEi2)
        lnTg2[:, i] .= t_axis

        # Pg
        y    .= t.lnPg[:, i]
        ip    = cubic_spline_interpolation(old_t, y, extrapolation_bc=Line())
        lnPg2[:, i] .= ip.(t_axis)

        # Ne
        y    .= t.lnNe[:, i]
        ip    = cubic_spline_interpolation(old_t, y, extrapolation_bc=Line())
        lnNe2[:, i] .= ip.(t_axis)

        # Ross
        y    .= t.lnRoss[:, i]
        ip    = cubic_spline_interpolation(old_t, y, extrapolation_bc=Line())
        lnRoss2[:, i] .= ip.(t_axis)
    end

    return RegularEoSTable(deepcopy(t.lnRho), lnTg2, lnEi2, lnPg2, lnRoss2, lnNe2)
end

lnT_to_lnEi(t) = error("No transformation from T to Ei implemented for this table.")

"""
Interpolate the function and first derivative of a cubic spline. (Taken from Tabgen)
    """
function interpolate_f_df!(ff, dd, xx, x, f, d)
    @inbounds for i in eachindex(ff)
        n = length(x)
        k = 0
        if (x[2]>x[1])
            for kin=2:n
                k = kin
                (xx[i] < x[kin]) && break
            end
        else
            for kin=2:n
                k = kin
                (x[kin] < xx[i]) && break
            end
        end

        begin 
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

"""
Add linear extrapolation outside of the given index.
"""
function add_outside!(f::AbstractArray{K, 1}, esi, eei) where {K}
    nF = length(f)
    for k=1:esi-1
        f[esi-k] = 2.0*f[esi] - f[esi+k]
    end 
    
    for k=1:nF-eei
        f[eei+k] = 2.0*f[eei] - f[eei-k]
    end 
end

"""
Compute the derivatives of all quantities w.r.t the axis variables. 
First: value, Second: df/dρ, third: df/dE (or dT)
"""
function derivative(t::RegularEoSTable)
    e_axis = ndims(t.lnT) == 1 ? false : true
    T = eltype(t.lnT)

    if e_axis
        lnDependent = zeros(T, size(t.lnPg)..., 3)
        lnAxis      = t.lnEi
        lnDependent[:, :, 1] .= t.lnT
    else
        lnDependent = zeros(T, size(t.lnPg)..., 3)
        lnAxis      = t.lnT
        lnDependent[:, :, 1] .= t.lnEi
    end

    nRhoBin  = length(t.lnRho)
    nAxisBin = length(lnAxis)

    lnPg = zeros(T, size(t.lnPg)..., 3)
    lnPg[:, :, 1] .= t.lnPg

    lnNe = zeros(T, size(t.lnPg)..., 3)
    lnNe[:, :, 1] .= t.lnNe

    lnRoss = zeros(T, size(t.lnPg)..., 3)
    lnRoss[:, :, 1] .= t.lnRoss

    for i=1:nAxisBin
        s = spline(interpolate(t.lnRho, lnDependent[i, :, 1], BSplineOrder(4)))
        lnDependent[i, :, 2] = (Derivative(1) * s).(t.lnRho)

        s = spline(interpolate(t.lnRho, lnPg[i, :, 1], BSplineOrder(4)))
        lnPg[i, :, 2] = (Derivative(1) * s).(t.lnRho)

        s = spline(interpolate(t.lnRho, lnNe[i, :, 1], BSplineOrder(4)))
        lnNe[i, :, 2] = (Derivative(1) * s).(t.lnRho)

        s = spline(interpolate(t.lnRho, lnRoss[i, :, 1], BSplineOrder(4)))
        lnRoss[i, :, 2] = (Derivative(1) * s).(t.lnRho)
    end

    for i=1:nRhoBin
        s   = spline(interpolate(lnAxis, lnDependent[:, i, 1], BSplineOrder(4)))
        lnDependent[:, i, 3] = (Derivative(1) * s).(lnAxis)

        s   = spline(interpolate(lnAxis, lnPg[:, i, 1], BSplineOrder(4)))
        lnPg[:, i, 3] = (Derivative(1) * s).(lnAxis)

        s   = spline(interpolate(lnAxis, lnNe[:, i, 1], BSplineOrder(4)))
        lnNe[:, i, 3] = (Derivative(1) * s).(lnAxis)

        s   = spline(interpolate(lnAxis, lnRoss[:, i, 1], BSplineOrder(4)))
        lnRoss[:, i, 3] = (Derivative(1) * s).(lnAxis)
    end 

    tab_new = if e_axis
        RegularEoSTable(t.lnRho, lnDependent, lnAxis, lnPg, lnRoss, lnNe)
    else
        RegularEoSTable(t.lnRho, lnAxis, lnDependent, lnPg, lnRoss, lnNe)
    end

    tab_new
end








#=========== Unify EoS and opacities on common internal Energy grid ==========#

function unify(eos::E, opacities::O, lnEi_new::AbstractVector{T2}) where {E<:EoSTable, O<:OpacityTable, T2<:AbstractFloat}
    nEiBin  = length(lnEi_new)
    nRhoBin = length(eos.lnRho)

    T = eltype(eos.lnRho)

    # New arrays
    lnTg2    = zeros(T, nEiBin, nRhoBin)
    lnNe2    = zeros(T, nEiBin, nRhoBin)
    lnRoss2  = zeros(T, nEiBin, nRhoBin)
    lnPg2    = zeros(T, nEiBin, nRhoBin)
    lnKross2 = zeros(T, nEiBin, nRhoBin)
    lnkappa2 = zeros(T, nEiBin, nRhoBin, length(opacities.λ))
    lnSrc2   = zeros(T, nEiBin, nRhoBin, length(opacities.λ))

    x = zeros(size(eos.lnEi, 1))
    y = zeros(size(eos.lnEi, 1))
    mask = zeros(Int, size(eos.lnEi, 1))

    @inbounds for i in 1:nRhoBin
        # This is a column of the table, interpolate it and evaluate at the right position

        # It is listed as a function of temperature -> get the new T axis for the new E axis
        mask .= sortperm(eos.lnEi[:, i])
        x    .= eos.lnEi[mask, i]
        y    .= eos.lnT[mask, i]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnTg2[:, i] .= ip.(lnEi_new)

        # Pg
        y    .= eos.lnPg[mask, i]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnPg2[:, i] .= ip.(lnEi_new)

        # Ne
        y    .= eos.lnNe[mask, i]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnNe2[:, i] .= ip.(lnEi_new)

        # Ross
        y    .= eos.lnRoss[mask, i]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnRoss2[:, i] .= ip.(lnEi_new)

        # kappa ross
        y    .= log.(opacities.κ_ross[mask, i])
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnKross2[:, i] .= ip.(lnEi_new)

        @inbounds for j in eachindex(opacities.λ)
            # kappa
            y    .= log.(opacities.κ[mask, i, j])
            ip    = linear_interpolation(x, y, extrapolation_bc=Line())
            lnkappa2[:, i, j] .= ip.(lnEi_new)

            # src
            y    .= log.(opacities.src[mask, i, j])
            ip    = linear_interpolation(x, y, extrapolation_bc=Line())
            lnSrc2[:, i, j] .= ip.(lnEi_new)
        end
    end

    return (RegularEoSTable(deepcopy(eos.lnRho), lnTg2, lnEi_new, lnPg2, lnRoss2, lnNe2), 
            RegularOpacityTable(exp.(lnkappa2), exp.(lnKross2),exp.(lnSrc2),deepcopy(opacities.λ), false))
end

function unify(eos_in::E, opacities_list::NTuple{N,O}; add_pad=0.0, upsample=-1) where {N, E<:EoSTable, O<:OpacityTable}
    eos = if ndims(eos_in.lnT) == 1
        @info "EoS passed that seems to be on the correct (lnT) grid. Puffing up..."
        puff_up(eos_in)
    else
        eos_in
    end

    aos = AxedEoS(eos_in)
    @assert DensityAxis(aos).dimension == 1

    emin = zeros(eltype(eos.lnRho), size(eos.lnEi, 2))
    emax = zeros(eltype(eos.lnRho), size(eos.lnEi, 2))
    for i in axes(eos.lnEi, 2)
        emin, emax = minimum(eos.lnEi[:, i]), maximum(eos.lnEi[:, i])
    end

    emin_pad = Base.convert(eltype(eos.lnRho), minimum(emin))
    emax_pad = Base.convert(eltype(eos.lnRho), maximum(emax))
    erange   = emax_pad - emin_pad
    emin_pad = Base.convert(eltype(eos.lnRho), emin_pad + abs(add_pad*erange))
    emax_pad = Base.convert(eltype(eos.lnRho), emax_pad - abs(add_pad*erange))

    new_axis = collect(range(emin_pad, emax_pad, length=upsample>0 ? upsample : size(eos.lnEi, 1)))
    
    @assert is_uniform(new_axis)

    opacities_out = []
    eos_out       = []
    for opacities in opacities_list
        e,o = unify(eos, opacities, new_axis)
        append!(opacities_out, [o])
        append!(eos_out, [e])
    end

    return if length(opacities_list)==1
        eos_out[1], opacities_out[1]
    else
        eos_out[1], (opacities_out...,)
    end
end

"""
    unify(eos, opacities; padding=0.001)

Interpolate EoS + opacities to a uniform internal energy grid. The tables must be on the same grid.
Nothing is assumed about the temperature grid. It only assumes that the density grid is structured.
Temperature can be unstructured. It will create the internal energy axis from the data and then will
interpolate as:

```jldoctest
for i in eachindex(eos.lnRho)
    ## Interpolate at constant rho
    x = eos.lnEi[:, i]
    y = eos.lnT[:, i]
    T_new[:, i] = linear_interpolation(x, y).(lnEi_new)

    ## and so on for the others...
end
```
"""
unify(eos::E, opa::O; kwargs...) where {E<:EoSTable, O<:OpacityTable} = unify(eos, (opa,); kwargs...)

"""
Combine the opacities from two opacity tables on a common EoS grid.
"""
function combine_opacities(eos1::E1, opacities1::O1, eos2::E2, opacities2::O2) where {E1<:EoSTable, O1<:OpacityTable, E2<:EoSTable, O2<:OpacityTable}
    eaxis = ndims(eos1.lnT)==1 ? false : true
    eaxis && @assert ndims(eos2.lnT)>1

    # decide on new grids
    lnVar_new = eaxis ? common_grid(eos1.lnEi,  eos2.lnEi) : common_grid(eos1.lnT,  eos2.lnT)
    lnRho_new = common_grid(eos1.lnRho, eos2.lnRho)

    do_regrid = if length(eos1.lnRho) == length(eos2.lnRho)
        if eaxis 
            (!all(eos1.lnEi .≈ eos2.lnEi)) | (!all(eos1.lnRho .≈ eos2.lnRho))
        else
            (!all(eos1.lnT .≈ eos2.lnT)) | (!all(eos1.lnRho .≈ eos2.lnRho))
        end
    else
        true
    end

    @info "Interpolation of tables to new, common grid needed: $(do_regrid)"

    # interpolate both tables to new grid 
    eos1_new, opacities1_new = do_regrid ? regrid(eos1, opacities1, lnRho_new, lnVar_new) : (eos1, opacities1)
    eos2_new, opacities2_new = do_regrid ? regrid(eos2, opacities2, lnRho_new, lnVar_new) : (eos2, opacities2)

    # Add the wavelength arrays together
    λ_new = vcat(opacities1_new.λ, opacities2_new.λ)
    umask = unique(i->λ_new[i], eachindex(λ_new))

    @info size(opacities1_new.λ) size(opacities2_new.λ)

    κ_new = zeros(eltype(opacities1_new.κ), size(opacities1_new.κ, 1), size(opacities1_new.κ, 2), length(λ_new)) #cat(opacities1_new.κ,   opacities2_new.κ,   dims=3)[:,:,umask]
    s_new = zeros(eltype(opacities1_new.κ), size(opacities1_new.κ, 1), size(opacities1_new.κ, 2), length(λ_new)) #cat(opacities1_new.src, opacities2_new.src, dims=3)[:,:,umask]

    #@inbounds for i in eachindex(opacities1_new.λ)
    #    κ_new[:, :, i] .= opacities1_new.κ[:, :, i]
    #    s_new[:, :, i] .= opacities1_new.src[:, :, i]
    #end

    #r = 1:length(opacities1_new.λ)
    ll = length(opacities1_new.λ) 
    κ_new[:, :, 1:ll] .= view(opacities1_new.κ,   :, :, 1:ll)
    s_new[:, :, 1:ll] .= view(opacities1_new.src, :, :, 1:ll)


    #ll = length(opacities1_new.λ)
    #@inbounds for i in eachindex(opacities2_new.λ)
    #    j = ll + i
    #    κ_new[:, :, j] .= opacities2_new.κ[:, :, i]
    #    s_new[:, :, j] .= opacities2_new.src[:, :, i]
    #end
   
    κ_new[:, :, (ll+1):(ll+length(opacities2_new.λ))] .= view(opacities2_new.κ,  :, :, 1:length(opacities2_new.λ))
    s_new[:, :, (ll+1):(ll+length(opacities2_new.λ))] .= view(opacities2_new.src,:, :, 1:length(opacities2_new.λ))

        
    κ_new   = κ_new[:, :, umask]
    GC.gc()
    s_new   = s_new[:, :, umask]
    GC.gc()
    λ_new   = λ_new[umask]
    GC.gc()

    # sort
    smask    = sortperm(λ_new)
    κ_new   .= view(κ_new, :, :, smask)
    s_new   .= view(s_new, :, :, smask)
    λ_new   .= view(λ_new, smask)

    combined_opacities = RegularOpacityTable(κ_new, deepcopy(opacities2_new.κ_ross), s_new, λ_new, opacities2_new.optical_depth)

    eos2_new, combined_opacities
end

"""
Get grid from 2 grids by picking the union
"""
common_grid(a1, a2) = begin    
    if length(a1) == length(a2)
        if all(a1 .≈ a2)
            return deepcopy(a1)
        end
    end

    aMax = min(a1[end], a2[end])
    aMin = min(a1[1],   a2[1])
    aLen = min(length(a1), length(a2))
    collect(range(aMin, aMax,length=aLen))
end


function cut(opacities::O, λ_lo, λ_hi) where {O<:OpacityTable}
    mask = (opacities.λ .>= λ_lo ) .& (opacities.λ .< λ_hi)
    RegularOpacityTable(opacities.κ[:, :, mask], deepcopy(opacities.κ_ross), opacities.src[:, :, mask], opacities.λ[mask], copy(opacities.optical_depth)) 
end







#================= Interaction of EoS with TS in/output ======================#

function write_as_stagger(lnT::Vector, lnRho::Vector; folder=@inWrapper("example/models"), teff=5777.0, logg=4.43, FeH=0.0) 
    sT = length(lnT)
    names_cols = String[]

    nbins = 1#ceil(Int, length(lnT) / 100)
    ibins = TSO.split_similar(lnT, nbins, mask=true)

    # write columns to file
    for column in eachindex(lnRho)
        for i in 1:nbins
            name = joinpath(folder, "TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")
            write_as_stagger(name; teff=teff, logg=logg, FeH=FeH, 
                                    T=exp.(lnT[ibins[i]]), rho=zeros(sT)[ibins[i]] .+ exp(lnRho[column]), id="TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")

            append!(names_cols, ["TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat"])
        end
    end

    open(joinpath(folder, "TSO_list.in"), "w") do f
        for name in names_cols
            write(f, name*"\n")
        end
    end
end

function write_as_stagger(lnT::Matrix, lnRho::Matrix; folder=@inWrapper("example/models"), teff=5777.0, logg=4.43, FeH=0.0) 
    sT = length(lnT)
    names_cols = String[]

    nbins = 1#ceil(Int, size(lnT, 1) / 100)
    ibins = TSO.split_similar(lnT[:, 1], nbins, mask=true)

    # write columns to file
    for column in axes(lnRho, 2)
        for i in 1:nbins
            name = joinpath(folder, "TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")
            write_as_stagger(name; teff=teff, logg=logg, FeH=FeH, 
                                    T=exp.(lnT[ibins[i], column]), rho=exp.(lnRho[ibins[i], column]), id="TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")

            append!(names_cols, ["TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat"])
        end
    end

    open(joinpath(folder, "TSO_list.in"), "w") do f
        for name in names_cols
            write(f, name*"\n")
        end
    end
end

"""
Read EoS columns from TS. Order: lnT, lnPe, lnRho, lnPg, lnEi, lnOp
"""
function read_eos_column(path, get_header=false)
    # Read the EoS column
    f = FortranFile(path, convert="big-endian")

    # Fist line with the size of the arrays
    # Order: lnT, lnPe, lnRho, lnPg, lnEi, lnOp
    ntau, xls = read(f, Int32, Float32)
    data = zeros(Float32, ntau, 6)

    for i in 1:ntau
        data[i, :] .= read(f, Float32, Float32, Float32, Float32, Float32, Float32)
    end

    close(f)

    return if get_header
        ntau, xls, data
    else
        data
    end
end

"""
Read Opacity columns from TS.
"""
function read_opacity_column(path)
    f = FortranFile(path, convert="big-endian")

    lheader = 0
    while true
        s = try
            Base.convert(String, read(f, FString{2}))
        catch
            " "
        end
        
        (s[end] !== '*') && break
        lheader += 1
    end

    # Skip the header
    FortranFiles.rewind(f)
    for i in 1:lheader
        read(f, FString{2})
    end

    abund = read(f, (Float32, 92))
    ntau, maxlam, λ_start, λ_end, Δλ = read(f, Int32, Int32, Float64, Float64, Float64)
    xls     = read(f, Float32)
    τ, vmic = read(f, (Float32, ntau), (Float32, ntau))
    κ_ross  = read(f, (Float32, ntau))
    κ_cont  = read(f, (Float32, ntau, maxlam))
    κ_line  = read(f, (Float32, ntau, maxlam))
    κ_scat  = read(f, (Float32, ntau, maxlam))
    close(f)

    opacity = zeros(Float32, ntau, maxlam, 3)
    opacity[:, :, 1] .= κ_cont
    opacity[:, :, 2] .= κ_line
    opacity[:, :, 3] .= κ_scat
    λ = Float32[λ_start:Δλ:λ_end...]

    (opacity, κ_ross, λ, abund)
end

matching_opacity(path) = begin
    new_p = path[first(findfirst("TSOeos", path)):findlast('.', path)-1]
    new_p = new_p[1:findlast('.', new_p)] *"dat.multi"
end

function _read_tables_eos_op(list_of_eos_tables, opacity_folder=@inTS(""); get_individual=false, kwargs...)
    # Get the dimensions
    nrows_per_table, ntables, T = get_dimensions(list_of_eos_tables)

    opacity_path = joinpath(opacity_folder, matching_opacity(list_of_eos_tables[1]))
  
    eaxis = false
    
    _,_,λ,abund = read_opacity_column(opacity_path)

    # Read the data
    lnρ    = zeros(T, ntables)
    lnT2   = zeros(T, nrows_per_table, ntables)
    lnPg   = zeros(T, nrows_per_table, ntables)
    lnEi2  = zeros(T, nrows_per_table, ntables)
    lnNe   = zeros(T, nrows_per_table, ntables)
    ross   = zeros(T, nrows_per_table, ntables)
    lnκ500 = zeros(T, nrows_per_table, ntables)
    κ      = zeros(T, nrows_per_table, ntables, length(λ))
    κc     = zeros(T, nrows_per_table, ntables, length(λ))
    κl     = zeros(T, nrows_per_table, ntables, length(λ))
    κs     = zeros(T, nrows_per_table, ntables, length(λ))
    Sν     = similar(κ)

    for (i,table) in enumerate(list_of_eos_tables)
        column, e_r = get_TSO_index(table)
        eos = read_eos_column(table)                                                   # The columns are ln(T), ln(pe), ln(ρ), ln(Pg), ln(E), ln(κ500)
        lnρ[column]        = eos[1, 3]
        lnT2[e_r, column]   .= eos[:, 1]
        lnPg[e_r, column]   .= eos[:, 4]
        lnEi2[e_r, column]  .= eos[:, 5]
        lnNe[e_r, column]   .= log.(exp.(eos[:, 2]) ./ (KBoltzmann .* exp.(eos[:, 1])))  # p = nkT
        lnκ500[e_r, column] .= eos[:, 6]                                                 # This is not the rosseland opacity! This is the 500nm for now

        opacity, κ_ross, _, _ = read_opacity_column(joinpath(opacity_folder, matching_opacity(table)))
        for j in axes(opacity, 1)
            opc = view(opacity, j, :, 1)
            opl = view(opacity, j, :, 2)
            ops = view(opacity, j, :, 3)

            opc[isnan.(opc) .| (opc .< 0.0)] .= 0.0
            opl[isnan.(opl) .| (opl .< 0.0)] .= 0.0
            ops[isnan.(ops) .| (ops .< 0.0)] .= 0.0
        end

        κ[e_r, column, :]  .= opacity[:, :, 1] .+ opacity[:, :, 2] .+ opacity[:, :, 3]
        κc[e_r, column, :] .= opacity[:, :, 1] 
        κl[e_r, column, :] .=                     opacity[:, :, 2] 
        κs[e_r, column, :] .=                                         opacity[:, :, 3]
        ross[e_r, column]  .= κ_ross
        Sν[e_r, column, :] .= Base.convert.(T, Bν(λ, exp.(lnT2[e_r, column]))) # Planck function (i.e. source function)
    end

    # check which of T/E is the axis and which is lnDependent
    for i in axes(lnT2, 2)
        if lnT2[:, i] != lnT2[:, 1]
            @info "Non-equal T columns detected. Assuming energy axis."
            eaxis = true
            break
        end
    end

    #=mask    = sortperm(lnρ)
    lnρ    .= lnρ[mask]
    lnPg   .= lnPg[:, mask]
    lnEi   .= lnEi[:, mask]
    lnNe   .= lnNe[:, mask]
    ross   .= ross[:, mask]
    κ      .= κ[:, mask, :] 
    lnκ500 .= lnκ500[:, mask]=#

    opacity_tables = if get_individual
        (RegularOpacityTable(κ, ross, Sν, λ, false), RegularOpacityTable(κc, ross, Sν, λ, false), RegularOpacityTable(κl, ross, Sν, λ, false), RegularOpacityTable(κs, ross, Sν, λ, false))
    else
        RegularOpacityTable(κ, ross, Sν, λ, false)
    end

    if eaxis
        (RegularEoSTable(lnρ, lnT2, lnEi2, lnPg, log.(ross), lnNe), opacity_tables)
    else
        (RegularEoSTable(lnρ, lnT2[:, 1], lnEi2, lnPg, log.(ross), lnNe), opacity_tables)
    end
end

function _read_tables_eos(list_of_eos_tables; kwargs...)
    # Get the dimensions
    nrows_per_table, ntables, T = get_dimensions(list_of_eos_tables)

    eaxis = false
    
    # Read the data
    lnρ    = zeros(T, ntables)
    lnT2   = zeros(T, nrows_per_table, ntables)
    lnPg   = zeros(T, nrows_per_table, ntables)
    lnEi2  = zeros(T, nrows_per_table, ntables)
    lnNe   = zeros(T, nrows_per_table, ntables)
    lnκ500 = zeros(T, nrows_per_table, ntables) 

    for (i,table) in enumerate(list_of_eos_tables)
        column, e_r = get_TSO_index(table)
        eos = read_eos_column(table)                                                   # The columns are ln(T), ln(pe), ln(ρ), ln(Pg), ln(E), ln(κ500)
        lnρ[column]          = eos[1, 3]
        lnT2[e_r, column]   .= eos[:, 1]
        lnPg[e_r, column]   .= eos[:, 4]
        lnEi2[e_r, column]  .= eos[:, 5]
        lnNe[e_r, column]   .= log.(exp.(eos[:, 2]) ./ (KBoltzmann .* exp.(eos[:, 1]))) # p = nkT
        lnκ500[e_r, column] .= eos[:, 6]                                                # This is not the rosseland opacity! This is the 500nm for now
    end

    # check which of T/E is the axis and which is lnDependent
    for i in axes(lnT2, 2)
        if lnT2[:, i] != lnT2[:, 1]
            @info "Non-equal T columns detected. Assuming energy axis."
            eaxis = true
            break
        end
    end

    #=mask    = sortperm(lnρ)
    lnρ    .= lnρ[mask]
    lnPg   .= lnPg[:, mask]
    lnEi   .= lnEi[:, mask]
    lnNe   .= lnNe[:, mask]
    ross   .= ross[:, mask]
    κ      .= κ[:, mask, :] 
    lnκ500 .= lnκ500[:, mask]=#

    if eaxis
        RegularEoSTable(lnρ, lnT2, lnEi2, lnPg, lnκ500, lnNe)
    else
        RegularEoSTable(lnρ, lnT2[:, 1], lnEi2, lnPg, lnκ500, lnNe)
    end
end

function get_dimensions(list_of_eos_tables)
    e = []
    d = []
    for i in eachindex(list_of_eos_tables)
        column, e_r = get_TSO_index(list_of_eos_tables[i])
        if column == 1
            append!(e, last(e_r))
        end

        append!(d, column)
    end

    eos = read_eos_column(list_of_eos_tables[1])
    #@show maximum(e), length(unique(d)), eltype(eos)
    return maximum(e), length(unique(d)), eltype(eos)
end

function read_tables(list_of_eos_tables::AbstractVector; kwargs...)
    opacity_path = @inTS "TSOeos_$(get_TSO_index(list_of_eos_tables[1])).multi"
    return if ispath(opacity_path)
        _read_tables_eos_op(list_of_eos_tables; kwargs...)
    else
        _read_tables_eos(list_of_eos_tables; kwargs...)
    end
end

read_tables(::Type{<:EoSTable}, ::Type{<:OpacityTable}, list_of_eos_tables::AbstractVector, opacity_folder=@inTS(""); kwargs...) = _read_tables_eos_op(list_of_eos_tables, opacity_folder; kwargs...)
read_tables(::Type{<:EoSTable}, list_of_eos_tables::AbstractVector; kwargs...) = _read_tables_eos(list_of_eos_tables; kwargs...)

read_tables(t::Type{<:EoSTable}, folder::String=".", args...; kwargs...) = read_tables(t, glob("_TSOeos_*_TSO.eos", folder), args...; kwargs...)
read_tables(t::Type{<:EoSTable}, t2::Type{<:OpacityTable}, folder::String=".", args...; kwargs...) = read_tables(t, t2, glob("_TSOeos_*_TSO.eos", folder), args...; kwargs...)

# User interface function
load(e::Type{<:EoSTable}, o::Type{<:OpacityTable}, args...; kwargs...) = read_tables(e, o, args...; kwargs...)
load(e::Type{<:EoSTable}, args...; kwargs...) = read_tables(e,    args...; kwargs...)
load(o::Type{<:OpacityTable}, e::Type{<:EoSTable}, args...; kwargs...) = begin
    e, o = read_tables(e, o, args...; kwargs...)
    o, e
end

get_TSO_index(name) = begin
    mask  = first(findfirst("TSOeos", name))
    maske = first(findlast(".", name)) -1
    rest  = name[mask:maske]
    rest  = rest[1:first(findlast(".", rest))-1]
    ident = split(split(rest, "_")[2], "-")
    d_col = parse(Int, ident[1])
    split(ident[2], ":")
    e_r   = range(parse(Int, split(ident[2], ":")[1]),parse(Int, split(ident[2], ":")[2]), step=1)

    d_col, e_r
end





#============================ Planck function ================================#

function Bλ(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0) #*aa_to_cm
    B[T .<= 1e-1] .= 0.0

    B
end

Bλ(λ::AbstractVector, T::AbstractVector) = begin
    B = zeros(length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= Bν(λ[i], T)
    end

    B
end

Bλ(λ::AbstractVector, T::AbstractMatrix) = begin
    B = zeros(size(T)..., length(λ))
    for i in axes(B, 3)
        B[:, :, i] .= Bν(λ[i], T)
    end

    B
end

function δBλ(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 * hc_k * exp(hc_k / (Λ*T)) / Λ^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2 #* aa_to_cm
    B[T .<= 1e-1] .= 0.0

    B
end

δBλ(λ::AbstractVector, T::AbstractVector) = begin
    B = zeros(length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= δBν(λ[i], T)
    end

    B
end

δBλ(λ::AbstractVector, T::AbstractMatrix) = begin
    B = zeros(size(T)..., length(λ))
    for i in axes(B, 3)
        B[:, :, i] .= δBν(λ[i], T)
    end

    B
end

@inline function δBλ(λ::A, T::A) where {A<:AbstractFloat}
    Λ = λ * aa_to_cm
    v = twohc2 * hc_k * exp(hc_k / (Λ*T)) / Λ^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2 #* aa_to_cm

    if isnan(v) || v<1e-30 
        A(1e-30)
    else
        v
    end
end

@inline function Bλ(λ::A, T::A) where {A<:AbstractFloat}
    Λ = λ * aa_to_cm
    v = twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0) #*aa_to_cm
    
    if isnan(v) || v<1e-30 
        A(1e-30)
    else
        v
    end
end

"""
Proxi for Bλ
"""
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





#=========================== Rosseland opacity ===============================#

"""
    new_rosseland_opacity(eos, opacities)

Integrate the rosseland opacity from the given monochromatic opacity table.
"""
function rosseland_opacity(eos::E, opacities; weights=ω_midpoint(opacities)) where {E<:AxedEoS}
    lnRoss = similar(eos.lnRoss)
    rosseland_opacity!(lnRoss, eos, opacities, weights=weights)
    lnRoss
end

"""
    new_rosseland_opacity(eos, opacities)

Integrate the rosseland opacity from the given monochromatic opacity table.
"""
function rosseland_opacity!(lnRoss, aos::E, opacities; weights=ω_midpoint(opacities)) where {E<:AxedEoS}
    eos      = aos.eos 
    axis_val = aos.energy_axes.values
    eaxis    = is_internal_energy(aos)
    
    lnRoss .= 0.0
    T       = zeros(eltype(eos.lnT), aos.energy_axes.length, aos.density_axes.length)
    if ndims(eos.lnT) == 2
        T .= exp.(eos.lnT)
    else
        puff_up!(T, exp.(eos.lnT))
    end

    B::Float32 = 0.0

    @inbounds for j in eachindex(eos.lnRho)
        @inbounds for i in eachindex(axis_val)
            B  = Float32(0.0)
            @inbounds for k in eachindex(opacities.λ)
                lnRoss[i, j] += weights[k] * 1 / opacities.κ[i, j, k] * δBν(opacities.λ[k], T[i, j])
                B += weights[k] * δBν(opacities.λ[k], T[i, j])
            end
            lnRoss[i, j] /= B
            lnRoss[i, j] = if 1.0 / lnRoss[i, j] == 0
                log(1e-30)
            else
                log(1.0 / lnRoss[i, j])
            end

            #if isnan(lnRoss[i, j])
            #    @warn "NaN in Ross - i:$(i), j:$(j), B:$(B), k:$(opacities.κ[i, j, 1])"
            #end
        end
    end
end

rosseland_opacity!(lnRoss, eos::E, args...; kwargs...) where {E<:RegularEoSTable} = rosseland_opacity!(lnRoss, AxedEoS(eos), args...; kwargs...) 
rosseland_opacity(eos::E, args...; kwargs...) where {E<:RegularEoSTable} = rosseland_opacity(AxedEoS(eos), args...; kwargs...)

"""
    transfer_rosseland(from, to)

Transfer rosseland opacity from "from" to "to". Can either be Opacity table or EoS.
"""
transfer_rosseland!(eos::E, opa::OpacityTable) where {E<:EoSTable} = opa.κ_ross .= exp.(eos.lnRoss);
transfer_rosseland!(opa::OpacityTable, eos::E) where {E<:EoSTable} = eos.lnRoss .= log.(opa.κ_ross);
transfer_rosseland!(eos::E, opa::OpacityTable) where {E<:AxedEoS}  = opa.κ_ross .= exp.(eos.eos.lnRoss);
transfer_rosseland!(opa::OpacityTable, eos::E) where {E<:AxedEoS}  = eos.eos.lnRoss .= log.(opa.κ_ross);






#===================== Optical depth related functions =======================#

"""
Compute monochromatic and rosseland optical depth scales of the model 
based on the opacity table
"""
function optical_depth(eos::E, opacities::OpacityTable, model::AbstractModel; binned=false) where {E<:AxedEoS}
    # Model z, ρ and T in cgs
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(opacities.κ)

    τ_λ    = zeros(T, length(lnρ), length(opacities.λ)) 
    τ_ross = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # For each wavelength we integrate along z, z[1]-> surface, z[end]-> bottom
    @inbounds for i in eachindex(opacities.λ)
        # Look up the opacity in the table
        κ  .= lookup(eos, opacities, :κ, lnρ, lnE, i)
        ρκ .= binned ? κ : exp.(lnρ) .* κ

        # Integrate: τ(z) = [ ∫ ρκ dz ]_z0 ^z
        @inbounds for j in eachindex(z)
            if j==1 
                τ_λ[1, i] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
            else
                τ_λ[j, i] = τ_λ[j-1, i] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
            end
        end
    end

    # Rosseland optical depth
    κ  .= lookup(eos, opacities, :κ_ross, lnρ, lnE)
    ρκ .= exp.(lnρ) .* κ # binned ? κ : exp.(lnρ) .* κ
    for j in eachindex(z)
        if j==1 
            τ_ross[1] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
        else
            τ_ross[j] = τ_ross[j-1] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
        end
    end

    return τ_ross, τ_λ
end

function rosseland_optical_depth(eos::E, opacities::OpacityTable, model::AbstractModel; binned=false) where {E<:AxedEoS}
    # Model z, ρ and T in cgs
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(opacities.κ)

    τ_ross = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # Rosseland optical depth
    κ  .= lookup(eos, opacities, :κ_ross, lnρ, lnE)
    ρκ .= exp.(lnρ) .* κ # binned ? κ : exp.(lnρ) .* κ
    for j in eachindex(z)
        if j==1 
            τ_ross[1] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
        else
            τ_ross[j] = τ_ross[j-1] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
        end
    end

    return τ_ross
end

function rosseland_optical_depth(eos::E, model::AbstractModel) where {E<:AxedEoS}
    # Model z, ρ and T in cgs
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos.eos.lnRoss)

    τ_ross = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # Rosseland optical depth
    κ  .= exp.(lookup(eos, :lnRoss, lnρ, lnE))
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

rosseland_optical_depth(eos::EoSTable, opacities::OpacityTable, model::AbstractModel; kwargs...)  = rosseland_optical_depth(@axed(eos), opacities, model; kwargs...)
rosseland_optical_depth(eos::EoSTable, model::AbstractModel; kwargs...)  = rosseland_optical_depth(@axed(eos), model; kwargs...)



function rosseland_depth(eos::E, opacities, model::AbstractModel) where {E<:AxedEoS}
    # Model z, ρ and T in cgs
    τ, lnρ, lnE = model.τ, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos.eos.lnRoss)

    z      = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # Rosseland optical depth
    κ  .= lookup(eos, opacities, :κ_ross, lnρ, lnE)
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

function rosseland_depth(eos::E, model::AbstractModel) where {E<:AxedEoS}
    # Model z, ρ and T in cgs
    τ, lnρ, lnE = model.τ, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(eos.eos.lnRoss)

    z      = zeros(T, length(lnρ)) 
    ρκ     = zeros(T, length(lnρ))
    κ      = zeros(T, length(lnρ))

    # Rosseland optical depth
    κ  .= exp.(lookup(eos, :lnRoss, lnρ, lnE))
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

rosseland_depth(eos::EoSTable, opacities::OpacityTable, model::AbstractModel; kwargs...)  = rosseland_depth(@axed(eos), opacities, model; kwargs...)
rosseland_depth(eos::EoSTable, model::AbstractModel; kwargs...)  = rosseland_depth(@axed(eos), model; kwargs...)







"""
Compute the formation height + opacity, i.e. the rosseland optical depth
where the monochromatic optical depth is 1.
"""
function formation_height(model::AbstractModel, eos::E, opacities::OpacityTable, τ_ross, τ_λ) where {E<:AxedEoS}
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT
    
    T = eltype(opacities.κ)

    z   = Base.convert.(T, z)
    lnρ = Base.convert.(T, lnρ)
    lnE = Base.convert.(T, lnE)

    rosseland_depth = zeros(T, size(τ_λ, 2))
    opacity_depth   = zeros(T, size(τ_λ, 2))

    lRoss = log.(τ_ross)
    lλ    = log.(τ_λ)

    t_mono = zeros(T, size(τ_λ, 1))
    r_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnρ, extrapolation_bc=Line())
    T_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnE, extrapolation_bc=Line())

    @inbounds for i in axes(τ_λ, 2)
        t_mono .= lλ[:, i]
        rosseland_depth[i] = exp.(linear_interpolation(Interpolations.deduplicate_knots!(t_mono, move_knots=true), lRoss, extrapolation_bc=Line())(0.0))
        opacity_depth[i] = lookup(eos, opacities, :κ, r_ross(rosseland_depth[i]), T_ross(rosseland_depth[i]), i)
    end

    rosseland_depth, opacity_depth
end

optical_depth(eos::RegularEoSTable,    args...; kwargs...) = optical_depth(AxedEoS(eos),    args...; kwargs...)
formation_height(model, eos::RegularEoSTable, args...; kwargs...) = formation_height(model, AxedEoS(eos), args...; kwargs...)

"""
Convert a monochromatic EoS into a binned EoS format
"""
toKappaLog!(opacities, eos) = begin
    opacities.src .= log.(opacities.src)

    for i in eachindex(eos.lnRho)
        opacities.κ[:, i, :] .*= exp(eos.lnRho[i])
    end
    opacities.κ .= log.(opacities.κ) 
end





#=================== Shift the Energy grid by a constant =====================#

"""
    shift_energy(eos, opacities, shift_amount)

Shift the energy grid by a constant. Because the spacing 
of the table has to be uniform in log, this requires 
    A) shifting the Energy grid as log(exp(e) + c)
    B) Interpolating to a new equidistant grid within the new limits.
"""
function shift_energy(aos::E, opacities, shift_amount) where {E<:AxedEoS}
    eos = aos.eos
    @assert aos.energy_axes == :lnEi

    eosNew = deepcopy(eos)
    eosNew.lnEi .= e_shifted.(eosNew.lnEi, shift_amount)

    new_Ei = range(first(eosNew.lnEi), last(eosNew.lnEi), length=length(eosNew.lnEi)) |> collect
    
    regrid(eosNew, opacities, eosNew.lnRho, new_Ei)
end

shift_energy(eos::RegularEoSTable, args...; kwargs...) = shift_energy(AxedEoS(eos), args...; kwargs...)

e_shifted(e, offset) = log(exp(e) + offset)








#============ replace the EoS corresponding to a given opacity table =========#

"""
    replace(eos_old, eos_new, opacities)

WARNING: The following is model dependent! It interpolates in E,
which may be different in different EoS. It should rather interpolate in T-rho

Replace an old EoS with a new one. Since the energy 
between the two might have been shifted, this function
proceeds as follows:
    A) use the old EoS to bisect the energy corresponding to
    the new eos lnT and lnRho
    B) Interpolate the old opacities on the old (regular) grid to the new grid in lnRho and lnEi
"""
function replaceEoS(eos_old::EoSTable, eos_new::EoSTable, opacities::OpacityTable)
    ## First we need to find the energy grid that corresponds to the new table energy "units"
    lnEi_newFromOld = similar(eos_new.lnT)
    emin,emax,_,_   = limits(eos_old)
    @inbounds for j in eachindex(eos_new.lnRho)
        @inbounds for i in eachindex(eos_new.lnEi)
            ## We search the Ei that will return, when called with rho, the temperature that is wanted
            lnEi_newFromOld[i, j] = bisect(eos_old, lnRho=eos_new.lnRho[j], lnT=eos_new.lnT[i, j], lnEi=[emin, emax])
        end
    end
    
    ## Now we know the the energy in the old units at every point in the new table
    ## So all that is left to do is to interpolate the old table and evaluate at 
    ## the new energy-density values
    κ      = similar(opacities.κ, size(eos_new.lnPg)..., length(opacities.λ))
    S      = similar(κ)
    κ_ross = similar(opacities.κ_ross, size(eos_new.lnPg)..., length(opacities.λ))
    
    #fr = table_interpolation(eos_old.lnEi, eos_old.lnRho, opacities.κ_ross)

    for k in eachindex(opacities.λ)
        fk = table_interpolation(eos_old.lnEi, eos_old.lnRho, view(opacities.κ,      :, :, k))
        fs = table_interpolation(eos_old.lnEi, eos_old.lnRho, view(opacities.src,    :, :, k))
        fr = table_interpolation(eos_old.lnEi, eos_old.lnRho, view(opacities.κ_ross, :, :, k))

        for j in eachindex(eos_new.lnRho)
            for i in eachindex(eos_new.lnEi)
                κ[i, j, k]      = fk(lnEi_newFromOld[i, j], eos_new.lnRho[j])
                S[i, j, k]      = fs(lnEi_newFromOld[i, j], eos_new.lnRho[j])
                κ_ross[i, j, k] = fr(lnEi_newFromOld[i, j], eos_new.lnRho[j])
            end
        end
    end

    #for j in eachindex(eos_new.lnRho)
    #    for i in eachindex(eos_new.lnEi)
    #        κ_ross[i, j] = fr(lnEi_newFromOld[i, j], os_new.lnRho[j])
    #    end
    #end

    RegularOpacityTable(κ, κ_ross, S, deepcopy(opacities.λ), opacities.optical_depth)
end 

@inline table_interpolation(x::A1, y::A2, z::A3) where {T<:AbstractFloat, A1<:AbstractArray{T, 1}, A2<:AbstractArray{T, 1}, A3<:AbstractArray{T, 2}} = begin
    extrapolate(interpolate((x, y), z, Gridded(Linear())), Line())                                                                                     
end





#====================== New Interpolation interface ==========================#

#= Please note that there is a new interpolation interface.                  =#
#= It does only partly rely on the lookup interface, however                 =#
#= it utilized the shape of the eos tables and should be much faster.        =#
#= It relys on dimensionwise 1D interpolation wherever possible.             =#

include("_interpolations.jl")

