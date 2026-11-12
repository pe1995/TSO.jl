# ============================================================================
# Functionality (no longer in use)
# ============================================================================

function transform(t::EoSTable, from_to::Pair{Symbol, Symbol})
    @warn "This method is no longer in use and should not be called."

    return if (first(from_to) == :lnT) & (last(from_to) == :lnEi)
        lnT_to_lnEi(t)
    else
        error("The given transformation of grids is not available yet.")
        nothing
    end
end

function lnT_to_lnEi(t::RegularEoSTable, interpolate_with_derivative=true)
    @warn "This method is no longer in use and should not be called."

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

function energy_grid(t::RegularEoSTable)
    @warn "This method is no longer in use and should not be called."

    T = eltype(t.lnEi)

    # new axis range
    lnEi2min = minimum(nanmin, t.lnEi)
    lnEi2max = maximum(nanmax, t.lnEi)

    # The upper limit should not correspond to a temperature larger than 500,000K (total upper limit)

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

function interpolate_f_df!(ff, dd, xx, x, f, d)
    @warn "This method is no longer in use and should not be called."

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

function add_outside!(f::AbstractArray{K, 1}, esi, eei) where {K}
    @warn "This method is no longer in use and should not be called."

    nF = length(f)
    for k=1:esi-1
        f[esi-k] = 2.0*f[esi] - f[esi+k]
    end

    for k=1:nF-eei
        f[eei+k] = 2.0*f[eei] - f[eei-k]
    end
end

function derivative(t::RegularEoSTable)
    @warn "This method is no longer in use and should not be called."

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

# ============================================================================
# Unify EoS and opacities on common internal Energy grid
# ============================================================================

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
        x    .= Interpolations.deduplicate_knots!(eos.lnEi[mask, i], move_knots=true)
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

unify(eos::E, opa::O; kwargs...) where {E<:EoSTable, O<:OpacityTable} = unify(eos, (opa,); kwargs...)

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

    ll = length(opacities1_new.λ)
    κ_new[:, :, 1:ll] .= view(opacities1_new.κ,   :, :, 1:ll)
    s_new[:, :, 1:ll] .= view(opacities1_new.src, :, :, 1:ll)

    #@inbounds for i in eachindex(opacities2_new.λ)
    #    κ_new[:, :, j] .= opacities2_new.κ[:, :, i]

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

# ============================================================================
# Interaction of EoS with TS in/output
# ============================================================================

include("_TS.jl")

# ============================================================================
# Planck function
# ============================================================================

function Bλ(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0) #*aa_to_cm
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
    B = @. twohc2 * hc_k * exp(hc_k / (Λ*T)) / Λ^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2 #* aa_to_cm
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
    new_rosseland_opacity(eos, opacities)

Integrate the rosseland opacity from the given monochromatic opacity table.
"""
rosseland_opacity(eos::E, opacities; weights=ω_midpoint(opacities)) where {E<:AxedEoS} = begin
    lnRoss = similar(eos.eos.lnRoss)
    rosseland_opacity!(lnRoss, eos, opacities, weights=weights)
    lnRoss
end

_rosseland_opacity!(lnRoss, lnRho, axis_val, l, weights, opacity, T) = _rosseland_opacityX!(lnRoss, lnRho, axis_val, l, weights, opacity, T)

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

_rosseland_opacity_chunk(chunk, lnRho, axis_val, T, l, weights, opacity) = begin
    BChunk = zeros(eltype(opacity), length(axis_val), length(lnRho))
    RossChunk = zeros(eltype(opacity), length(axis_val), length(lnRho))
    dbnuChunk = Ref(eltype(opacity)(0.0))

    @inbounds for k in chunk
        @inbounds for j in eachindex(lnRho)
            @inbounds for i in eachindex(axis_val)
                δBν!(dbnuChunk, l[k], T[i, j])
                val = weights[k] * dbnuChunk[]
                if !isnan(val) && val > 0.0
                    op = opacity[i, j, k]
                    if op > 0.0
                        RossChunk[i, j] += val / op
                        BChunk[i, j] += val
                    end
                end
            end
        end
    end

    return [BChunk, RossChunk]
end

_rosseland_opacity2!(lnRoss, lnRho, axis_val, l, weights, opacity, db) = begin
    B::Float32 = Float32(0.0)

    @inbounds for j in eachindex(lnRho)
        @inbounds for i in eachindex(axis_val)
            B = Float32(0.0)

            lnRoss[i, j] = sum(weights .* 1 ./ view(opacity, :, i, j) .* view(db, :, i, j))
            B = sum(weights .* view(db, :, i, j))

            lnRoss[i, j] /= B
            lnRoss[i, j] = if 1.0 / lnRoss[i, j] == 0
                NaN
            else
                log(1.0 / lnRoss[i, j])
            end
        end
    end
end

rosseland_opacity!(lnRoss, eos::E, args...; kwargs...) where {E<:RegularEoSTable} = rosseland_opacity!(lnRoss, AxedEoS(eos), args...; kwargs...)
rosseland_opacity(eos::E, args...; kwargs...) where {E<:RegularEoSTable} = rosseland_opacity(AxedEoS(eos), args...; kwargs...)

"""
    rosseland_opacity!(lnRoss, aos::E, opacities; weights=ω_midpoint(opacities))

Integrate the rosseland opacity from the given monochromatic opacity table.
"""
function rosseland_opacity!(lnRoss, aos::E, opacities; weights=ω_midpoint(opacities)) where {E<:AxedEoS}
    eos      = aos.eos
    axis_val = aos.energy_axes.values
    eaxis    = is_internal_energy(aos)

    T = zeros(eltype(eos.lnT), aos.energy_axes.length, aos.density_axes.length)
    if ndims(eos.lnT) == 2
        T .= exp.(eos.lnT)
    else
        puff_up!(T, exp.(eos.lnT))
    end

    @optionalTiming rosseland_time _rosseland_opacityX!(
        lnRoss, eos.lnRho, axis_val, opacities.λ, weights, opacities.κ, T
    )

end

# ============================================================================
# transfer Rosseland opacity between EoS and OpacityTable
# ============================================================================

"""
    transfer_rosseland(from, to)

Transfer rosseland opacity from "from" to "to". Can either be Opacity table or EoS.
"""
transfer_rosseland!(eos::E, opa::OpacityTable) where {E<:EoSTable} = opa.κ_ross .= exp.(eos.lnRoss);
transfer_rosseland!(opa::OpacityTable, eos::E) where {E<:EoSTable} = eos.lnRoss .= log.(opa.κ_ross);
transfer_rosseland!(eos::E, opa::OpacityTable) where {E<:AxedEoS}  = opa.κ_ross .= exp.(eos.eos.lnRoss);
transfer_rosseland!(opa::OpacityTable, eos::E) where {E<:AxedEoS}  = eos.eos.lnRoss .= log.(opa.κ_ross);

# ============================================================================
# Optical depth related functions
# ============================================================================

_optical_depth!(ρκ, τ, z, λ, lnρ, lnE, eos, opacities, binned) = begin
    eosx, opax = extended(eos), extended(opacities)
    cRho, cE = weights(eosx, lnρ, lnE; extrapolate=true)
    # For each wavelength we integrate along z, z[1]-> surface, z[end]-> bottom
    @inbounds for i in eachindex(λ)
        sample!((ρκ,), eosx, opax, (:κ,), cRho, cE, i)
        ρκ .= binned ? ρκ : exp.(lnρ) .* ρκ

        # Integrate: τ(z) = [ ∫ ρκ dz ]_z0 ^z
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

_optical_depthX!(τ, z, λ, lnρ, lnE, eos, opacities, binned) = begin
    chunks = Iterators.partition(eachindex(λ), length(λ) ÷ Threads.nthreads()) |> collect

    tasks = map(chunks) do chunk
        Threads.@spawn _optical_depth_chunk!(chunk, z, λ, lnρ, lnE, eos, opacities, binned)
    end
    results = fetch.(tasks)

    for (i, chunk) in enumerate(chunks)
        τ[:, chunk] .= results[i]
    end
end

_optical_depth_chunk!(chunk, z, λ, lnρ, lnE, eos, opacities, binned) = begin
    τ = zeros(eltype(opacities.κ), length(z), length(chunk))
    ρκ = zeros(eltype(opacities.κ), length(z))
    ρ = exp.(lnρ)
    eosx, opax = extended(eos), extended(opacities)
    cRho, cE = weights(eosx, lnρ, lnE; extrapolate=true)

    # For each wavelength we integrate along z, z[1]-> surface, z[end]-> bottom
    @inbounds for (i, k) in enumerate(chunk)
        sample!((ρκ,), eosx, opax, (:κ,), cRho, cE, k)
        ρκ .= binned ? ρκ : ρ .* ρκ

        # Integrate: τ(z) = [ ∫ ρκ dz ]_z0 ^z
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
    optical_depth(eos, opacities, model; binned=false)

Compute monochromatic and rosseland optical depth scales of the model
based on the opacity table.
"""
function optical_depth(eos::E, opacities::OpacityTable, model::AbstractModel; binned=false) where {E<:AxedEoS}
    # make sure model is in right direction
    model = TSO.flip(model, depth=true)

    # Model z, ρ and T in cgs
    z, lnρ, lnE = model.z, model.lnρ, is_internal_energy(eos) ? model.lnEi : model.lnT

    T = eltype(opacities.κ)

    τ_λ = zeros(T, length(lnρ), length(opacities.λ))
    τ_ross = zeros(T, length(lnρ))

    # For each wavelength we integrate along z, z[1]-> surface, z[end]-> bottom
    @optionalTiming optical_depth_time if Threads.nthreads() <= 1
        ρκ = zeros(T, length(lnρ))
        _optical_depth!(ρκ, τ_λ, z, opacities.λ, lnρ, lnE, eos, opacities, binned)
    else
        _optical_depthX!(τ_λ, z, opacities.λ, lnρ, lnE, eos, opacities, binned)
    end

    # Rosseland optical depth
    τ_ross .= rosseland_optical_depth(eos, opacities, model, binned=binned)

    # impose limits
    τ_λ[τ_λ .> 1e30] .= 1e30
    τ_λ[τ_λ .< 1e-30] .= 1e-30

    τ_ross[τ_ross .> 1e30] .= 1e30
    τ_ross[τ_ross .< 1e-30] .= 1e-30

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
    κ  .= sample(eos, opacities, (:κ_ross,), lnρ, lnE; extrapolate=true)[1]
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
    κ  .= exp.(sample(eos, (:lnRoss,), lnρ, lnE; extrapolate=true)[1])
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
    κ  .= sample(eos, opacities, (:κ_ross,), lnρ, lnE; extrapolate=true)[1]
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
    κ  .= exp.(sample(eos, (:lnRoss,), lnρ, lnE; extrapolate=true)[1])
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

# ============================================================================
# formation height
# ============================================================================

"""
    formation_height(model, eos, opacities, τ_ross, τ_λ)

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

    lRoss = log10.(τ_ross)
    lλ    = log10.(τ_λ)

    λ = opacities.λ

    r_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnρ, extrapolation_bc=Flat())
    T_ross = linear_interpolation(Interpolations.deduplicate_knots!(lRoss, move_knots=true), lnE, extrapolation_bc=Flat())

    if Threads.nthreads() > 1
        @info "Computing formation height with $(Threads.nthreads()) threads."
    end

    @optionalTiming formation_height_time begin
        chunks = Iterators.partition(eachindex(λ), length(λ) ÷ Threads.nthreads()) |> collect

        tasks = map(chunks) do chunk
            Threads.@spawn _formation_height_chunk(chunk, τ_λ, lλ, lRoss, r_ross, T_ross, eos, opacities)
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

_formation_height_chunk(chunk, τ_λ, lλ, lRoss, r_ross, T_ross, eos, opacities) = begin
    T = eltype(opacities.κ)
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

        opacity_depth[i] = sample(eos, opacities, (:κ,), r_ross(rosseland_depth[i]), T_ross(rosseland_depth[i]), k; extrapolate=true)[1]
    end

    return [rosseland_depth, opacity_depth]
end

optical_depth(eos::RegularEoSTable,    args...; kwargs...) = optical_depth(AxedEoS(eos),    args...; kwargs...)
formation_height(model, eos::RegularEoSTable, args...; kwargs...) = formation_height(model, AxedEoS(eos), args...; kwargs...)

toKappaLog!(opacities, eos) = begin
    opacities.src .= log.(opacities.src)

    for i in eachindex(eos.lnRho)
        opacities.κ[:, i, :] .*= exp(eos.lnRho[i])
    end
    opacities.κ .= log.(opacities.κ)
end

# ============================================================================
# Shift the Energy grid by a constant
# ============================================================================

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

# ============================================================================
# replace the EoS corresponding to a given opacity table
# ============================================================================

function replaceEoS(eos_old::EoSTable, eos_new::EoSTable, opacities::OpacityTable)
    ## First we need to find the energy grid that corresponds to the new table energy "units"
    lnEi_newFromOld = similar(eos_new.lnT)
    emin,emax,_,_   = limits(eos_old)
    @inbounds for j in eachindex(eos_new.lnRho)
        @inbounds for i in eachindex(eos_new.lnEi)
            ## We search the Ei that will return, when called with rho, the temperature that is wanted
            lnEi_newFromOld[i, j] = sample(eos_old, :lnEi; lnRho=eos_new.lnRho[j], lnT=eos_new.lnT[i, j])
        end
    end

    ## Now we know the the energy in the old units at every point in the new table
    ## So all that is left to do is to interpolate the old table and evaluate at
    ## the new energy-density values
    κ      = similar(opacities.κ, size(eos_new.lnPg)..., length(opacities.λ))
    S      = similar(κ)
    κ_ross = similar(opacities.κ_ross, size(eos_new.lnPg)..., length(opacities.λ))

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

    #        κ_ross[i, j] = fr(lnEi_newFromOld[i, j], os_new.lnRho[j])

    RegularOpacityTable(κ, κ_ross, S, deepcopy(opacities.λ), opacities.optical_depth)
end

@inline table_interpolation(x::A1, y::A2, z::A3) where {T<:AbstractFloat, A1<:AbstractArray{T, 1}, A2<:AbstractArray{T, 1}, A3<:AbstractArray{T, 2}} = begin
    extrapolate(interpolate((x, y), z, Gridded(Linear())), Line())
end

# ============================================================================
# New Interpolation interface
# Please note that there is a new interpolation interface.
# It does only partly rely on the lookup interface, however
# it utilized the shape of the eos tables and should be much faster.
# It relys on dimensionwise 1D interpolation wherever possible.
# ============================================================================

include("_interpolations.jl")

# ============================================================================
# Timer setup
# ============================================================================
const rosseland_time = Ref(false)
const optical_depth_time = Ref(false)
const formation_height_time = Ref(false)
