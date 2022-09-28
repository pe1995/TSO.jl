
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

    @show lnEi2min lnEi2max

    lnT   = t.lnT
    lnRho = t.lnRho

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
        lneo .= t.lnEi[:, i, 1]         # ln ei on lnT grid
        f    .= lnT                     # lnT grid values
        d    .= 1.0 ./ t.lnEi[:, i, 3]  # d lnT / d ln(ei) | rho on lnT grid
        g    .= t.lnEi[:, i, 2]         # d ln ei / d ln(rho)| Tg on lnT grid
       
        # get h = d/d ln ei ( d ln ei / d ln rho) on lnT grid
        @show i
        if i == 16
            @show lneo g typeof(lneo)
        end
        s   = spline(interpolate(lneo, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lneo)
        
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
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

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
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

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

        # the rosseland opacity ###############################################
        f .= t.lnRoss[:, i, 1] # f on lnT grid
        d .= t.lnRoss[:, i, 3] # d f / d ln tg on lnT grid
        g .= t.lnRoss[:, i, 2] # d f / d ln rho on lnT grid

        # get h= d/d ln ei ( d f / d ln rho) on lnT grid
        s   = spline(interpolate(lnT, g, BSplineOrder(4)))
        h   = (Derivative(1) * s).(lnT)
        
        interpolate_f_df!(view(ff, esi:eei), view(dd, esi:eei), lnT2, lnT, f, d)
        interpolate_f_df!(view(gg, esi:eei), view(hh, esi:eei), lnT2, lnT, g, h)

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

lnT_to_lnEi(t) = error("No transformation from T to Ei implemented for this table.")

"""Interpolate the function and first derivative of a cubic spline. (Taken from Tabgen)"""
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

"""Add linear extrapolation outside of the given index."""
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
    
    if e_axis
        lnDependent = zeros(size(t.lnPg)..., 3)
        lnAxis      = t.lnEi
        lnDependent[:, :, 1] .= t.lnT
    else
        lnDependent = zeros(size(t.lnPg)..., 3)
        lnAxis      = t.lnT
        lnDependent[:, :, 1] .= t.lnEi
    end

    nRhoBin  = length(t.lnRho)
    nAxisBin = length(lnAxis)

    lnPg = zeros(size(t.lnPg)..., 3)
    lnPg[:, :, 1] .= t.lnPg

    lnNe = zeros(size(t.lnPg)..., 3)
    lnNe[:, :, 1] .= t.lnNe

    lnRoss = zeros(size(t.lnPg)..., 3)
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

function write_as_stagger(lnT::Vector, lnRho::Vector; folder=@inWrapper("example/models"), teff=5777.0, logg=4.43, FeH=0.0) 
    sT = length(lnT)
    names_cols = String[]

    # write columns to file
    for column in eachindex(lnRho)
        name = joinpath(folder, "TSOeos_$(column).dat")
        write_as_stagger(name; teff=teff, logg=logg, FeH=FeH, 
                                T=exp.(lnT), rho=zeros(sT) .+ exp(lnRho[column]), id="TSOeos_$(column)")

        append!(names_cols, ["TSOeos_$(column).dat"])
    end

    open(joinpath(folder, "TSO_list.in"), "w") do f
        for name in names_cols
            write(f, name*"\n")
        end
    end
end

"""Read EoS columns from TS. Order: lnT, lnPe, lnRho, lnPg, lnEi, lnOp"""
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

"""Read Opacity columns from TS."""
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

function read_tables(list_of_eos_tables)
    # Get the dimensions
    eos_col = read_eos_column(list_of_eos_tables[1])
    T = eltype(eos_col)

    nrows_per_table = size(eos_col, 1)
    ntables         = length(list_of_eos_tables)
    @assert ntables == nrows_per_table

    _,_,λ,abund = read_opacity_column(@inTS("TSOeos_$(get_TSO_index(list_of_eos_tables[1])).multi"))

    lnT = eos_col[:, 1]
    B   = Base.convert.(T, Bν(λ, exp.(lnT))) # Planck function (i.e. source function)

    # Read the data
    lnρ    = zeros(T, ntables)
    lnPg   = zeros(T, nrows_per_table, ntables)
    lnEi   = zeros(T, nrows_per_table, ntables)
    lnNe   = zeros(T, nrows_per_table, ntables)
    ross   = zeros(T, nrows_per_table, ntables)
    lnκ500 = zeros(T, nrows_per_table, ntables)
    κ      = zeros(T, nrows_per_table, ntables, length(λ))
    Sν     = similar(κ)

    for (i,table) in enumerate(list_of_eos_tables)
        eos = read_eos_column(table)          # The columns are ln(T), ln(pe), ln(ρ), ln(Pg), ln(E), ln(κ500)
        lnρ[i]        = eos[1, 3]
        lnPg[:, i]   .= eos[:, 4]
        lnEi[:, i]   .= eos[:, 5]
        lnNe[:, i]   .= log.(exp.(eos[:, 2]) ./ (1.38065e-16 .* exp.(eos[:, 1]))) # p = nkT
        lnκ500[:, i] .= eos[:, 6] # This is not the rosseland opacity! This is the 500nm for now

        opacity, κ_ross, _, _ = read_opacity_column(@inTS("TSOeos_$(get_TSO_index(table)).multi"))
        for j in 1:nrows_per_table
            opc = view(opacity, j, :, 1)
            opl = view(opacity, j, :, 2)
            ops = view(opacity, j, :, 3)

            opc[isnan.(opc)] .= 0.0
            opl[isnan.(opl)] .= 0.0
            ops[isnan.(ops)] .= 0.0
        end

        κ[:, i, :]  .= opacity[:, :, 1] .+ opacity[:, :, 2] .+ opacity[:, :, 3]
        ross[:, i]  .= κ_ross
        Sν[:, i, :] .= B
    end

    mask    = sortperm(lnρ)
    lnρ    .= lnρ[mask]
    lnPg   .= lnPg[:, mask]
    lnEi   .= lnEi[:, mask]
    lnNe   .= lnNe[:, mask]
    ross   .= ross[:, mask]
    κ      .= κ[:, mask, :] 
    lnκ500 .= lnκ500[:, mask]

    (RegularEoSTable(lnρ, lnT, lnEi, lnPg, log.(ross), lnNe), RegularOpacityTable(κ, ross, Sν, λ, false))
end

get_TSO_index(name) = begin
    mask = first(findfirst("TSOeos", name))
    rest = name[mask:end]
    parse(Int, split(rest, "_")[2])
end

function Bν(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 /Λ^5 /(exp(hc_k / (Λ*T)) - 1.0) *aa_to_cm
    B[T .<= 1e-1] .= 0.0

    B
end

Bν(λ::AbstractArray, T::AbstractArray) = begin
    B = zeros(length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= Bν(λ[i], T)
    end

    B
end

function δBν(λ::AbstractFloat, T::AbstractArray)
    Λ = λ * aa_to_cm
    B = @. twohc2 * hc_k * exp(hc_k / (Λ*T)) / la^6 / T^2 / (exp(hc_k / (Λ*T))-1)^2 * AA_TO_CM
    B[T .<= 1e-1] .= 0.0

    B
end

δBν(λ::AbstractArray, T::AbstractArray) = begin
    B = zeros(length(T), length(λ))
    for i in axes(B, 2)
        B[:, i] .= δBν(λ[i], T)
    end

    B
end