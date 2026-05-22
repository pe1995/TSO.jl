#= Interpolate tables =#

is_scattered(ax::EoSAxis) = ax.dimension > 1
is_gridded(ax::EoSAxis)   = !is_scattered(ax)
is_gridded(aos::AxedEoS)  = is_gridded(DensityAxis(aos)) & is_gridded(EnergyAxis(aos))


#= 
    If at least one of the axes is gridded, things can go faster,
    so the first method should always be "make this gridded, without changing the axes".
=#


"""
    gridded(eos, opacities...)

Interpolate a random eos, opacity table to a grid. Will depending on the input
decide what interpolation is needed to save time. If switch=true, additionally the energy axis will be switched.
"""
function gridded(aos::A, opacities::OpacityTable...; conservative=false, upsample=-1, switch=false, newE=nothing, newRho=nothing) where {A<:AxedEoS}
    # Check the current grid, if it is gridded alredy there is nothing to do
    if is_gridded(aos) & !switch ## Nothing to do here
        return (aos.eos, opacities...)  
    #elseif switch
    #    return switch_energy(aos, opacities...; conservative=conservative, upsample=upsample, newE=nothing)
    end

    # check if at least parts are gridded
    density_gridded = is_gridded(DensityAxis(aos))
    energy_gridded  = is_gridded(EnergyAxis(aos))

    newE2, newRho2 = pick_axis(aos, conservative=conservative, upsample=upsample, switch=switch)
    newE   = isnothing(newE)   ? newE2   : newE
    newRho = isnothing(newRho) ? newRho2 : newRho

    @assert is_uniform(newE) 
    @assert is_uniform(newRho)

    ename = switch ? dependent_energy_variable(aos) : energy_variable(aos)

    # Pick the best method for interpolation, worst case pick the generic one
    eos_new, opacities_new = if density_gridded      # Only the density is gridded already 
        interpolate_at_density(aos, opacities...; (ename=>newE,)...)
    elseif energy_gridded                            # Only the energy is gridded already
        interpolate_at_energy(aos, opacities...; lnRho=newRho)
    else                                             # Nothing is gridded
        interpolate_2D(aos, opacities...; lnRho=newRho, (ename=>newE,)...)
    end

    (eos_new, opacities_new...)
end

"""
    uniform(eos, opacities...)

Make the EoS uniform in both Axes. Checks if it is gridded before.
"""
function uniform(aos::A, opacities::OpacityTable...; conservative=false, upsample=-1, switch=false, newE=nothing, newRho=nothing) where {A<:AxedEoS}
    if is_uniform(aos) & !switch ## Nothing to do here
        return (aos.eos, opacities...)  
    #elseif switch
    #    return switch_energy(aos, opacities...; conservative=conservative, upsample=upsample, newE=nothing)
    end

    density_gridded = is_gridded(DensityAxis(aos))
    energy_gridded  = is_gridded(EnergyAxis(aos))
    density_uniform = is_uniform(DensityAxis(aos))
    energy_uniform  = is_uniform(EnergyAxis(aos))

    newE2, newRho2 = pick_axis(aos, conservative=conservative, upsample=upsample, switch=switch)
    newE   = isnothing(newE)   ? newE2   : newE
    newRho = isnothing(newRho) ? newRho2 : newRho

    @assert is_uniform(newE) 
    @assert is_uniform(newRho)

    ename = switch ? dependent_energy_variable(aos) : energy_variable(aos)

    @show density_gridded
    @show energy_gridded 
    @show density_uniform
    @show energy_uniform 
    @show size(newE)

    return if density_uniform                  ## This mean E is not uniform!
        aos_puff = if energy_gridded           ## However E is still gridded, which means we have to puff it up in order
            puff_up(aos, :energy)
        else
            aos
        end
        interpolate_at_density(aos_puff, opacities...; (ename=>newE,)...)

    elseif energy_uniform
        aos_puff = if density_gridded          ## However d is still gridded, which means we have to puff it up in order
            puff_up(aos, :density)
        else
            aos
        end
        interpolate_at_energy(aos_puff, opacities...; lnRho=newRho)

    else ## They are both not uniform, so only the order of interpolation needs to be decided 
        if energy_gridded
            aos_puff = if density_gridded          ## However d is still gridded, which means we have to puff it up in order
                puff_up(aos, :density)
            else
                aos
            end
            aos, opacities... = interpolate_at_energy(aos_puff, opacities...; lnRho=newRho)

            aos_puff = puff_up(AxedEoS(aos), :energy)
            interpolate_at_density(aos_puff, opacities...; (ename=>newE,)...)
        elseif density_gridded
            aos_puff = if energy_gridded           ## However E is still gridded, which means we have to puff it up in order
                puff_up(aos, :energy)
            else
                aos
            end
            aos, opacities... = interpolate_at_density(aos_puff, opacities...; (ename=>newE,)...)

            aos_puff = puff_up(AxedEoS(aos), :density)
            interpolate_at_energy(aos_puff, opacities...; lnRho=newRho)
        else
            interpolate_2D(aos, opacities...; lnRho=newRho, (ename=>newE,)...)
        end
    end
end    

"""
    complement(eos_old, eos_new, opacities...)

Switch out the EoS from all opacities with the new one.
"""
function complement(aos_old::E, aos_new::E2, opacities::OpacityTable...) where {E<:AxedEoS, E2<:AxedEoS} 
    @assert is_gridded(aos_new)

    ename_new = EnergyAxis(aos_new).name
    evals_new = EnergyAxis(aos_new).values

    density_gridded = is_gridded(DensityAxis(aos_old))
    energy_gridded  = is_gridded(EnergyAxis(aos_old))

    ## we only need to make sure that if the new axis is internal energy we dont have any offsets
    aos_old_mod = if is_internal_energy(aos_new)
        ee, rr      = meshgrid(EnergyAxis(aos_old).values, DensityAxis(aos_old).values)
        told        = is_internal_energy(aos_old) ? aos_old.eos.lnT : ee
        newEOldGrid = lookup(aos_new, :lnEi, rr, told)
        
        rho_mod = if is_internal_energy(aos_old)
            rr
        else
            deepcopy(aos_old.eos.lnRho)
        end

        e_old_mod = RegularEoSTable(rho_mod, aos_old.eos.lnT, newEOldGrid, 
                                    aos_old.eos.lnPg, aos_old.eos.lnRoss, aos_old.eos.lnNe)
        
        eaxis = if is_internal_energy(aos_old) # we have replaced the axis
            EnergyAxis(e_old_mod, axis=:lnEi)
        else
            aos_old.energy_axes
        end 

        daxis = if is_internal_energy(aos_old) # we have replaced the axis
            DensityAxis(e_old_mod)
        else
            aos_old.density_axes
        end

        AxedEoS(e_old_mod, eaxis, daxis)
    else
        aos_old
    end

    e = aos_new.eos

    density_gridded = is_gridded(DensityAxis(aos_old_mod))
    energy_gridded  = is_gridded(EnergyAxis(aos_old_mod))

    # Interpolate the old tables to the grid
    efinal, ofinal... = interpolate_2D(aos_old_mod, opacities...; lnRho=deepcopy(e.lnRho), (ename_new=>evals_new,)...)
    return if length(ofinal) == 1
        first(ofinal)
    else
        ofinal
    end
end

"""
    switch_energy(eos, opacities; kwargs...)

Switch out the Energy axis of the given tables. It is required that the tables are gridded.
"""
function switch_energy(aos::A, opacities::OpacityTable...; conservative=false, upsample=-1, switch=true, newE=nothing, newRho=nothing) where {A<:AxedEoS} 
    @assert is_gridded(aos)

    newE2,_ = pick_axis(aos, conservative=conservative, upsample=upsample, switch=true)
    newE = isnothing(newE) ? newE2 : newE

    ename = dependent_energy_variable(aos)

    ## this is easy because we know it is gridded
    aos_puff = puff_up(aos, :energy)
    interpolate_at_density(aos_puff, opacities...; (ename=>newE,)...)
end

gridded(eos::A, args...; kwargs...) where {A<:RegularEoSTable} = gridded(AxedEoS(eos), args...; kwargs...)
uniform(eos::A, args...; kwargs...) where {A<:RegularEoSTable} = uniform(AxedEoS(eos), args...; kwargs...)
complement(eos_old::E, eos_new::E2, opacities::OpacityTable...) where {E<:RegularEoSTable, E2<:RegularEoSTable} = complement(AxedEoS(eos_old), AxedEoS(eos_new), opacities...)
switch_energy(eos::A, args...; kwargs...) where {A<:RegularEoSTable} = switch_energy(AxedEoS(eos), args...; kwargs...)





#================================================================= Interpolation functions ===#

"""
    interpolate_at_density(eos, opacities...; newE...)

Interpolate the eos + opacities assuming a grid in density. 
This is eqivalent, but more general version of the old "unify" function.
This function is only called if the Energy axis is NOT gridded, i.e. it is always 2D.
"""
function interpolate_at_density(aos::E, opacities...; newE...) where {E<:AxedEoS}
    @assert !is_gridded(EnergyAxis(aos)) & is_gridded(DensityAxis(aos))

    eos = aos.eos
    newAxis_name = first(keys(newE))
    newAxis_val  = newE[newAxis_name]

    oldVname = newAxis_name == :lnEi ? :lnT : :lnEi
    
    nEBin   = length(newAxis_val)
    nRhoBin = length(eos.lnRho)

    T = eltype(eos.lnRho)

    # New arrays
    lnV2     = zeros(T, nEBin, nRhoBin)
    lnNe2    = zeros(T, nEBin, nRhoBin)
    lnRoss2  = zeros(T, nEBin, nRhoBin)
    lnPg2    = zeros(T, nEBin, nRhoBin)

    refo     = first(opacities)
    lnKross2 = zeros(T, nEBin, nRhoBin, length(opacities))
    lnkappa2 = zeros(T, nEBin, nRhoBin, length(refo.λ), length(opacities))
    lnSrc2   = zeros(T, nEBin, nRhoBin, length(refo.λ), length(opacities))

    x = zeros(size(eos.lnPg, 1))
    y = zeros(size(eos.lnPg, 1))
    mask = zeros(Int, size(eos.lnPg, 1))

    newAxisOldE = getfield(aos.eos, newAxis_name)
    oldlnV      = getfield(aos.eos, oldVname)

    @inbounds for i in 1:nRhoBin
        # This is a column of the table, interpolate it and evaluate at the right position
        mask .= sortperm(view(newAxisOldE, :, i))
        x    .= Interpolations.deduplicate_knots!(newAxisOldE[mask, i], move_knots=true)
        y    .= @view(oldlnV[mask, i])
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnV2[:, i] .= ip.(newAxis_val)

        # Pg
        y    .= @view(eos.lnPg[mask, i])
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnPg2[:, i] .= ip.(newAxis_val)

        # Ne
        y    .= @view(eos.lnNe[mask, i])
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnNe2[:, i] .= ip.(newAxis_val)

        # Ross
        y    .= @view(eos.lnRoss[mask, i])
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnRoss2[:, i] .= ip.(newAxis_val)

        # kappa ross
        for m in eachindex(opacities)
            opa   = opacities[m]
            y    .= log.(@view(opa.κ_ross[mask, i]))
            ip    = linear_interpolation(x, y, extrapolation_bc=Line())
            lnKross2[:, i, m] .= ip.(newAxis_val)

            @inbounds for j in eachindex(opa.λ)
                # kappa
                y    .= log.(@view(opa.κ[mask, i, j]))
                ip    = linear_interpolation(x, y, extrapolation_bc=Line())
                lnkappa2[:, i, j, m] .= ip.(newAxis_val)

                # src
                y    .= log.(@view(opa.src[mask, i, j]))
                ip    = linear_interpolation(x, y, extrapolation_bc=Line())
                lnSrc2[:, i, j, m] .= ip.(newAxis_val)
            end
        end
    end

    return if newAxis_name == :lnEi
        (RegularEoSTable(deepcopy(eos.lnRho), lnV2, newAxis_val, lnPg2, lnRoss2, lnNe2), 
        (RegularOpacityTable(exp.(lnkappa2[:, :, :, i]), exp.(lnKross2[:, :, i]), exp.(lnSrc2[:, :, :, i]), deepcopy(opacities[i].λ), opacities[i].optical_depth) for i in eachindex(opacities))...)
    else
        (RegularEoSTable(deepcopy(eos.lnRho), newAxis_val, lnV2, lnPg2, lnRoss2, lnNe2), 
        (RegularOpacityTable(exp.(lnkappa2[:, :, :, i]), exp.(lnKross2[:, :, i]), exp.(lnSrc2[:, :, :, i]), deepcopy(opacities[i].λ), opacities[i].optical_depth) for i in eachindex(opacities))...)
    end
end

"""
    interpolate_at_energy(eos, opacities...; lnRho)

Interpolate the eos + opacities assuming a grid in energy. 
This function is only called if the Density axis is NOT gridded, i.e. it is always 2D.
"""
function interpolate_at_energy(aos::E, opacities...; lnRho) where {E<:AxedEoS}
    @assert is_gridded(EnergyAxis(aos)) & !is_gridded(DensityAxis(aos))

    eos = aos.eos
       
    nEBin   = EnergyAxis(aos).length
    nRhoBin = length(lnRho)

    T = eltype(eos.lnPg)

    # New arrays
    lnV2     = zeros(T, nEBin, nRhoBin)
    lnNe2    = zeros(T, nEBin, nRhoBin)
    lnRoss2  = zeros(T, nEBin, nRhoBin)
    lnPg2    = zeros(T, nEBin, nRhoBin)

    refo     = first(opacities)
    lnKross2 = zeros(T, nEBin, nRhoBin, length(opacities))
    lnkappa2 = zeros(T, nEBin, nRhoBin, length(refo.λ), length(opacities))
    lnSrc2   = zeros(T, nEBin, nRhoBin, length(refo.λ), length(opacities))

    x = zeros(size(eos.lnPg, 2))
    y = zeros(size(eos.lnPg, 2))
    mask = zeros(Int, size(eos.lnPg, 2))

    r_old  = eos.lnRho
    oldlnV = getfield(eos, dependent_energy_variable(aos)) # We need to lookup only the thing which is not the axis

    @inbounds for i in 1:nEBin
        # This is a column of the table, interpolate it and evaluate at the right position
        mask .= sortperm(view(r_old, i, :))
        x    .= Interpolations.deduplicate_knots!(r_old[i, mask])
        y    .= @view oldlnV[i, mask]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnV2[i, :] .= ip.(lnRho)

        # Pg
        y    .= @view eos.lnPg[i, mask]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnPg2[i, :] .= ip.(lnRho)

        # Ne
        y    .= @view eos.lnNe[i, mask]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnNe2[i, :] .= ip.(lnRho)

        # Ross
        y    .= @view eos.lnRoss[i, mask]
        ip    = linear_interpolation(x, y, extrapolation_bc=Line())
        lnRoss2[i, :] .= ip.(lnRho)

        # kappa ross
        for m in eachindex(opacities)
            opa   = opacities[m]
            y    .= log.(@view(opa.κ_ross[i, mask]))
            ip    = linear_interpolation(x, y, extrapolation_bc=Line())
            lnKross2[i, :, m] .= ip.(lnRho)

            @inbounds for j in eachindex(opa.λ)
                # kappa
                y    .= log.(@view(opa.κ[i, mask, j]))
                ip    = linear_interpolation(x, y, extrapolation_bc=Line())
                lnkappa2[i, :, j, m] .= ip.(lnRho)

                # src
                y    .= log.(@view(opa.src[i, mask, j]))
                ip    = linear_interpolation(x, y, extrapolation_bc=Line())
                lnSrc2[i, :, j, m] .= ip.(lnRho)
            end
        end
    end

    return if energy_variable(aos) == :lnEi
        (RegularEoSTable(lnRho, lnV2, deepcopy(EnergyAxis(aos).values), lnPg2, lnRoss2, lnNe2), 
        (RegularOpacityTable(exp.(lnkappa2[:, :, :, i]), exp.(lnKross2[:, :, i]), exp.(lnSrc2[:, :, :, i]), deepcopy(opacities[i].λ), opacities[i].optical_depth) for i in eachindex(opacities))...)
    else
        (RegularEoSTable(lnRho, deepcopy(EnergyAxis(aos).values), lnV2, lnPg2, lnRoss2, lnNe2), 
        (RegularOpacityTable(exp.(lnkappa2[:, :, :, i]), exp.(lnKross2[:, :, i]), exp.(lnSrc2[:, :, :, i]), deepcopy(opacities[i].λ), opacities[i].optical_depth) for i in eachindex(opacities))...)
    end
end

"""
    interpolate_2D(eos, opacities...; new_axis...)

Interpolate the eos + opacities assuming nothing is gridded.
TODO: When extrapolation happens this is will return NaN values. There
needs to be interpolation in both directions, first Rho and E, where
NaN values are filled.
"""
function interpolate_2D(aos::E, opacities...; lnRho, newE...) where {E<:AxedEoS}
    skip_opacities = length(opacities) == 0
    eos = aos.eos
    newAxis_name = first(keys(newE))
    newAxis      = newE[newAxis_name]
    newRho       = lnRho

    @assert is_uniform(newAxis)
    @assert is_uniform(newRho)

    eos = aos.eos
    iseaxis  = is_internal_energy(aos)
    var_name = dependent_energy_variable(aos)

    newVar  = zeros(eltype(eos.lnPg), length(newAxis), length(newRho))
    newPg   = zeros(eltype(eos.lnPg), length(newAxis), length(newRho))
    newNe   = zeros(eltype(eos.lnPg), length(newAxis), length(newRho))
    newRoss = zeros(eltype(eos.lnPg), length(newAxis), length(newRho))

    PgLook      = lookup_function(aos, :lnPg)
    NeLook      = lookup_function(aos, :lnNe)
    RossLook    = lookup_function(aos, :lnRoss)
    VarLook     = lookup_function(aos, var_name)
    oldAxisLook = lookup_function(aos, aos.energy_axes.name)

    old_is_new = newAxis_name == energy_variable(aos)

    # if the axis has been switched, we can not lookup all variables at the new grid points
    # so in this case we need to lookup the actual energy variable first for the new point
    # and use this to lookup the rest, otherwise we can lookup everythinf directly
    oldAxisNewGrid = zeros(eltype(eos.lnPg), length(newAxis), length(newRho))

    @inbounds for j in eachindex(newRho)
        @inbounds for i in eachindex(newAxis)
            # if the new axis is the same as the old axis, we can just lookup the new axis directly
            oldAxisNewGrid[i, j] = old_is_new ? newAxis[i] : lookup(oldAxisLook, newRho[j], newAxis[i])
            newPg[i, j]          = lookup(PgLook,      newRho[j], oldAxisNewGrid[i, j])
            newNe[i, j]          = lookup(NeLook,      newRho[j], oldAxisNewGrid[i, j])
            newRoss[i, j]        = lookup(RossLook,    newRho[j], oldAxisNewGrid[i, j])
            newVar[i, j]         = lookup(VarLook,     newRho[j], oldAxisNewGrid[i, j])
        end
    end

    # For this it might actually be smarter to use a meshgrid
    κ_all, s_all, κ_rall = if !skip_opacities
        ee, rr = meshgrid(newAxis, newRho)

        ref_o  = first(opacities)
        κ_all  = zeros(eltype(ref_o.κ), length(newAxis), length(newRho), length(ref_o.λ), length(opacities))
        s_all  = zeros(eltype(ref_o.κ), length(newAxis), length(newRho), length(ref_o.λ), length(opacities))
        κ_rall = zeros(eltype(ref_o.κ), length(newAxis), length(newRho), length(opacities))

        kLF = [lookup_function(aos, opacities[m], :κ_ross) for m in eachindex(opacities)]

        for m in eachindex(opacities)
            κ_rall[:, :, m] = lookup(kLF[m], rr, oldAxisNewGrid)
        end

        for k in eachindex(ref_o.λ)
            kLF = [lookup_function(aos, opacities[m], :κ, k)   for m in eachindex(opacities)]
            sLF = [lookup_function(aos, opacities[m], :src, k) for m in eachindex(opacities)]

            for m in eachindex(opacities)
                κ_all[:, :, k, m] .= lookup(kLF[m], rr, oldAxisNewGrid)
                s_all[:, :, k, m] .= lookup(sLF[m], rr, oldAxisNewGrid)
            end
        end

        κ_all, s_all, κ_rall 
    else
        nothing, nothing, nothing
    end

    # old E - new T -> oldAxisNewGrid - E -> newVar - T
    # old E - new E -> oldAxisNewGrid - E -> newVar - T
    # old T - new T -> oldAxisNewGrid - T -> newVar - E
    # old T - new E -> oldAxisNewGrid - T -> newVar - E
    newT = if iseaxis & old_is_new # old E - new E
        newVar
    elseif iseaxis & !old_is_new   # old E - new T
        newAxis
    elseif !iseaxis & old_is_new   # old T - new T
        newAxis
    else                           # old T - new E
        oldAxisNewGrid
    end

    newE = if iseaxis & old_is_new # old E - new E
        newAxis
    elseif iseaxis & !old_is_new   # old E - new T
        oldAxisNewGrid
    elseif !iseaxis & old_is_new   # old T - new T
        newVar
    else                           # old T - new E
        newAxis
    end

    return if !skip_opacities
        eos, opacities... = (RegularEoSTable(newRho, newT, newE, newPg, newRoss, newNe), (RegularOpacityTable(κ_all[:, :, :, m], κ_rall[:, :, m], s_all[:, :, :, m], opacities[m].λ, opacities[m].optical_depth) for m in eachindex(opacities))...)
        fill_nan!(AxedEoS(eos), opacities...)

        (eos, opacities...)
    else
        eos = RegularEoSTable(newRho, newT, newE, newPg, newRoss, newNe)
        fill_nan!(eos)

        eos
    end
end





#============================================================ EoS manipulation functions ===#

function upsample(eos, axis, N)
    a = getfield(eos, axis)
    new_axis = range(minimum(a), maximum(a), length=N) |> collect

    replace_axis(eos; (axis=>new_axis,)...)
end

"""
    replace_axis(eos; new_axis=new_axis)

Replace a axis of the given EoS with the new one.
Name of the new axis has to match one of the current axes of the table.
"""
function replace_axis(aos::E; axis...) where {E<:AxedEoS}
    eos   = aos.eos
    eaxis = aos.energy_axes
    daxis = aos.density_axes

    new_size = [eaxis.length, daxis.length]
    for axis_name in keys(axis)
        if !((axis_name == eaxis.name) | (axis_name == daxis.name))
            error("Given new axis is not an axis of the given EoS.")
        end

        if axis_name == eaxis.name
            new_size[1] = length(axis[axis_name])
        elseif axis_name == daxis.name
            new_size[2] = length(axis[axis_name])
        end
    end

    lnRho_new = zeros(eltype(aos.eos.lnPg), new_size[2])
    lnAxi_new = zeros(eltype(aos.eos.lnPg), new_size[1])

    for axis_name in keys(axis)
        if axis_name == eaxis.name
            lnAxi_new .= axis[axis_name]
        elseif axis_name == daxis.name
            lnRho_new .= axis[axis_name]
        end
    end

    if all(i->i==0, lnRho_new)
        lnRho_new .= daxis.values
    end
    if all(i->i==0, lnAxi_new)
        lnAxi_new .= eaxis.values
    end

    iseaxis  = is_internal_energy(aos)
    var_name = iseaxis ? :lnT : :lnEi

    newVar  = zeros(eltype(eos.lnRho), length(lnAxi_new), length(lnRho_new))
    newPg   = zeros(eltype(eos.lnRho), length(lnAxi_new), length(lnRho_new))
    newNe   = zeros(eltype(eos.lnRho), length(lnAxi_new), length(lnRho_new))
    newRoss = zeros(eltype(eos.lnRho), length(lnAxi_new), length(lnRho_new))

    PgLook   = lookup_function(eos, :lnPg)
    NeLook   = lookup_function(eos, :lnNe)
    RossLook = lookup_function(eos, :lnRoss)
    VarLook  = lookup_function(eos, var_name)

    for j in eachindex(lnRho_new)
        for i in eachindex(lnAxi_new)
            newPg[i, j]   = lookup(PgLook,   lnRho_new[j], lnAxi_new[i])
            newNe[i, j]   = lookup(NeLook,   lnRho_new[j], lnAxi_new[i])
            newRoss[i, j] = lookup(RossLook, lnRho_new[j], lnAxi_new[i])
            newVar[i, j]  = lookup(VarLook,  lnRho_new[j], lnAxi_new[i])
        end
    end

    iseaxis ? RegularEoSTable(lnRho_new, newVar, lnAxi_new, newPg, newRoss, newNe) :
              RegularEoSTable(lnRho_new, lnAxi_new, newVar, newPg, newRoss, newNe)
end

replace_axis(eos::E; kwargs...) where {E<:RegularEoSTable} = replace_axis(AxedEoS(eos); kwargs...)





#================================================================== Utilities ===#

conservative_axis(mat::Vector) = default_axis(mat)
conservative_axis(mat::Matrix) = begin
    mirow = zeros(eltype(mat), size(mat, 1))
    marow = zeros(eltype(mat), size(mat, 1))

    for j in axes(mat, 1)
        mirow[j] = minimum(view(mat, j, :))
        marow[j] = maximum(view(mat, j, :))
    end

    micol = zeros(eltype(mat), size(mat, 2))
    macol = zeros(eltype(mat), size(mat, 2))

    for j in axes(mat, 2)
        micol[j] = minimum(view(mat, :, j))
        macol[j] = maximum(view(mat, :, j))
    end

    sort([max(maximum(mirow),maximum(micol)), min(minimum(marow), minimum(macol))])
end

default_axis(mat) = (minimum(mat), maximum(mat))

pick_axis(eos; conservative=false, upsample=-1, switch=true) = begin
    switch_values = if switch
        is_internal_energy(eos) ? getfield(eos.eos, :lnT) : getfield(eos.eos, :lnEi)
    else
        EnergyAxis(eos).values
    end

    newAxis = if conservative
        upsample==-1 ? range(conservative_axis(switch_values)..., length=EnergyAxis(eos).length) :
                       range(conservative_axis(switch_values)..., length=upsample)
    else    
        upsample==-1 ? range(default_axis(switch_values)..., length=EnergyAxis(eos).length) :
                       range(default_axis(switch_values)..., length=upsample)
    end

    newRho = if conservative & !is_uniform(DensityAxis(eos))
        range(conservative_axis(DensityAxis(eos).values)..., length=DensityAxis(eos).length) 
    elseif !is_uniform(DensityAxis(eos))
        range(default_axis(DensityAxis(eos).values)..., length=DensityAxis(eos).length) 
    else
        deepcopy(DensityAxis(eos).values)
    end

    sort(collect(newAxis)), sort(collect(newRho))
end


"""
    puff_up(eos)

Puff up the temperature axis to a 2D array.
"""
function puff_up(eos::RegularEoSTable)
    if ndims(eos.lnT) == 2
        return deepcopy(eos)
    end

    T = similar(eos.lnPg)
    puff_up!(T, eos.lnT, 1)

    RegularEoSTable(deepcopy(eos.lnRho), T, deepcopy(eos.lnEi), deepcopy(eos.lnPg), deepcopy(eos.lnRoss), deepcopy(eos.lnNe))
end

"""
    puff_up!(Tnew, Told)

Puff up the temperature axis to a 2D array.
"""
function puff_up!(Tnew::AbstractArray{T,2}, Told::Vector, dim=1) where {T}
    @inbounds for j in axes(Tnew, 2)
        @inbounds for i in axes(Tnew, 1)
            Tnew[i, j] = dim==1 ? Told[i] : Told[j]
        end
    end
end

puff_up(aos::AxedEoS, which::Symbol) = begin
    return if which == :energy
        puff_up_energy(aos)
    elseif which == :density
        puff_up_density(aos)
    else
        error("Pick energy or density.")
    end
end

puff_up_energy(aos::AxedEoS) = begin
    E_puff = similar(aos.eos.lnPg)     ## to use the functions
    puff_up!(E_puff, EnergyAxis(aos).values, 1) 

    newEi = is_internal_energy(aos) ? E_puff      : deepcopy(aos.eos.lnEi)
    newT  = is_internal_energy(aos) ? deepcopy(aos.eos.lnT) : E_puff

    @axed RegularEoSTable(deepcopy(aos.eos.lnRho), newT, newEi, deepcopy(aos.eos.lnPg), deepcopy(aos.eos.lnRoss), deepcopy(aos.eos.lnNe))
end

puff_up_density(aos::AxedEoS) = begin
    d_puff = similar(aos.eos.lnPg)     ## to use the functions
    puff_up!(d_puff, DensityAxis(aos).values, 2) 
    @axed RegularEoSTable(d_puff, deepcopy(aos.eos.lnT), deepcopy(aos.eos.lnEi), deepcopy(aos.eos.lnPg), deepcopy(aos.eos.lnRoss), deepcopy(aos.eos.lnNe))
end




#=================================================================== NaN handling ===#

@inline nan_or_inf(a) = isnan(a) | !isfinite(a)
@inline interpolate_at(x, y, x0; bc=Flat()) = linear_interpolation(Interpolations.deduplicate_knots!(deepcopy(x), move_knots=true), y, extrapolation_bc=bc).(x0)
"""
    _lerp_flat(x_knots, y_knots, n, x_query)

Zero-allocation piecewise linear interpolation with flat (constant) extrapolation.
`n` is the number of valid knots in the buffer arrays.
"""
@inline function _lerp_flat(x_knots, y_knots, n, x_query)
    @inbounds x_query <= x_knots[1] && return y_knots[1]
    @inbounds x_query >= x_knots[n] && return y_knots[n]

    # Binary search for the interval [lo, lo+1] containing x_query
    lo, hi = 1, n
    @inbounds while hi - lo > 1
        mid = (lo + hi) >> 1
        if x_knots[mid] <= x_query
            lo = mid
        else
            hi = mid
        end
    end

    @inbounds begin
        t = (x_query - x_knots[lo]) / (x_knots[hi] - x_knots[lo])
        return y_knots[lo] + t * (y_knots[hi] - y_knots[lo])
    end
end

"""
    _interpolate_at!(y_slice, x_coords, row_mask, buf_x, buf_y; is_log, fallback)

Interpolate NaN/Inf values in `y_slice` using the coordinate array `x_coords`.
Operates in-place with pre-allocated buffers (`buf_x`, `buf_y`) to avoid allocations.
If `is_log=true`, the data is interpolated in log space (input/output stays linear).
Returns `true` if any NaN/Inf was found (and `row_mask` is filled accordingly).
"""
function _interpolate_at!(y_slice, x_coords, row_mask, buf_x, buf_y;
                           is_log=false, fallback=1.0f-30)
    # Detect NaN/Inf (and non-positive values when working in log space)
    @inbounds for i in eachindex(y_slice)
        if is_log
            row_mask[i] = nan_or_inf(y_slice[i]) || y_slice[i] <= 0
        else
            row_mask[i] = nan_or_inf(y_slice[i])
        end
    end
    n_bad = count(row_mask)
    n_bad == 0 && return false

    n_good = length(row_mask) - n_bad
    if n_good >= 2
        # Pack good points into pre-allocated buffers
        k = 0
        @inbounds for i in eachindex(x_coords)
            if !row_mask[i]
                k += 1
                buf_x[k] = x_coords[i]
                buf_y[k] = is_log ? log(y_slice[i]) : y_slice[i]
            end
        end

        @inbounds for i in eachindex(x_coords)
            if row_mask[i]
                val = _lerp_flat(buf_x, buf_y, k, x_coords[i])
                y_slice[i] = is_log ? exp(val) : val
            end
        end
    else
        # Too few good points — use fallback
        @inbounds for i in eachindex(y_slice)
            if row_mask[i]
                y_slice[i] = fallback
            end
        end
    end
    return true
end

"""
    fill_nan!(aos, opacities...; fallback=1.0f-30)

Repair NaN/Inf values in the EoS and opacity tables by linear interpolation.
Two passes are performed:
  1. At constant ρ (column), interpolate across the energy/T axis — primary direction.
  2. At constant T (row), interpolate across log ρ — fallback for points that
     could not be repaired in pass 1 (fewer than 2 good neighbours in T).

Returns a `Matrix{Bool}` of size `(nE, nρ)` that is `true` at every `(T, ρ)` grid
point where *any* field originally contained NaN or Inf.
"""
function fill_nan!(aos::A, opacities::OpacityTable...; fallback=1.0f-30) where {A<:AxedEoS}
    eos   = aos.eos
    x_rho = DensityAxis(aos).values     
    x_e   = EnergyAxis(aos).values      
    nE    = length(x_e)
    nρ    = length(x_rho)

    nan_mask = falses(nE, nρ)           
    var      = getfield(eos, dependent_energy_variable(aos))
    lfallback = log(fallback)

    has_opacities = length(opacities) > 0

    # Pass 1: constant ρ, interpolate across T/energy axis
    col_mask = Vector{Bool}(undef, nE)
    buf_xe   = similar(x_e)
    buf_ye   = similar(x_e)

    @inbounds for j in 1:nρ
        for field in (view(var, :, j), view(eos.lnPg, :, j),
                      view(eos.lnNe, :, j), view(eos.lnRoss, :, j))
            had = _interpolate_at!(field, x_e, col_mask, buf_xe, buf_ye;
                                    is_log=false, fallback=lfallback)
            had && (view(nan_mask, :, j) .|= col_mask)
        end

        if has_opacities
            for opa in opacities
                had = _interpolate_at!(view(opa.κ_ross, :, j), x_e, col_mask,
                                       buf_xe, buf_ye; is_log=true, fallback=fallback)
                had && (view(nan_mask, :, j) .|= col_mask)

                @inbounds for k in eachindex(opa.λ)
                    had = _interpolate_at!(view(opa.κ, :, j, k), x_e, col_mask,
                                           buf_xe, buf_ye; is_log=true, fallback=fallback)
                    had && (view(nan_mask, :, j) .|= col_mask)

                    if size(opa.src) != (1, 1, 1)
                        had = _interpolate_at!(view(opa.src, :, j, k), x_e, col_mask,
                                               buf_xe, buf_ye; is_log=true, fallback=fallback)
                        had && (view(nan_mask, :, j) .|= col_mask)
                    end
                end
            end
        end
    end

    # Pass 2: constant T, interpolate across ρ (fallback)
    row_mask = Vector{Bool}(undef, nρ)
    buf_x    = similar(x_rho)
    buf_y    = similar(x_rho)

    @inbounds for i in 1:nE
        for field in (view(var, i, :), view(eos.lnPg, i, :),
                      view(eos.lnNe, i, :), view(eos.lnRoss, i, :))
            _interpolate_at!(field, x_rho, row_mask, buf_x, buf_y;
                              is_log=false, fallback=lfallback)
        end

        if has_opacities
            for opa in opacities
                _interpolate_at!(view(opa.κ_ross, i, :), x_rho, row_mask,
                                  buf_x, buf_y; is_log=true, fallback=fallback)

                @inbounds for k in eachindex(opa.λ)
                    _interpolate_at!(view(opa.κ, i, :, k), x_rho, row_mask,
                                      buf_x, buf_y; is_log=true, fallback=fallback)

                    if size(opa.src) != (1, 1, 1)
                        _interpolate_at!(view(opa.src, i, :, k), x_rho, row_mask,
                                          buf_x, buf_y; is_log=true, fallback=fallback)
                    end
                end
            end
        end
    end

    return Base.convert(Array{Bool, 2}, nan_mask)
end

fill_nan!(eos::RegularEoSTable, args...; kwargs...) = fill_nan!(AxedEoS(eos), args...; kwargs...)





function set_small!(aos, opa, small=1e-30)
    otherV = view(getfield(aos.eos, dependent_energy_variable(aos)), :, :)
    lsmall = log(small)
    @inbounds for j in eachindex(aos.density_axes.values)
        @inbounds for i in eachindex(aos.energy_axes.values)
            otherV[i, j]         = check(otherV[i, j], lsmall) ? otherV[i, j] : lsmall
            aos.eos.lnPg[i, j]   = check(aos.eos.lnPg[i, j],   lsmall) ? aos.eos.lnPg[i, j]   : lsmall
            aos.eos.lnNe[i, j]   = check(aos.eos.lnNe[i, j],   lsmall) ? aos.eos.lnNe[i, j]   : lsmall
            aos.eos.lnRoss[i, j] = check(aos.eos.lnRoss[i, j], lsmall) ? aos.eos.lnRoss[i, j] : lsmall

            if !isnothing(opa)
                opa.κ_ross[i, j]     = check(opa.κ_ross[i, j],     small)  ? opa.κ_ross[i, j]     : small
                has_src = size(opa.src) != (1, 1, 1)
                @inbounds for k in eachindex(opa.λ)
                    opa.κ[i, j, k]     = check(opa.κ[i, j, k],   small)  ? opa.κ[i, j, k]    : small
                    if has_src
                        opa.src[i, j, k]   = check(opa.src[i, j, k], small)  ? opa.src[i, j, k]  : small
                    end
                end
            end
        end
    end
end

function set_large!(aos, opa, large=1e30)
    otherV = view(getfield(aos.eos, dependent_energy_variable(aos)), :, :)
    llarge = log(large)
    @inbounds for j in eachindex(aos.density_axes.values)
        @inbounds for i in eachindex(aos.energy_axes.values)
            otherV[i, j]         = check_large(otherV[i, j], llarge) ? otherV[i, j] : llarge
            aos.eos.lnPg[i, j]   = check_large(aos.eos.lnPg[i, j],   llarge) ? aos.eos.lnPg[i, j]   : llarge
            aos.eos.lnNe[i, j]   = check_large(aos.eos.lnNe[i, j],   llarge) ? aos.eos.lnNe[i, j]   : llarge
            aos.eos.lnRoss[i, j] = check_large(aos.eos.lnRoss[i, j], llarge) ? aos.eos.lnRoss[i, j] : llarge
            
            if !isnothing(opa)
                opa.κ_ross[i, j]     = check_large(opa.κ_ross[i, j],     large)  ? opa.κ_ross[i, j]     : large
                has_src = size(opa.src) != (1, 1, 1)
                @inbounds for k in eachindex(opa.λ)
                    opa.κ[i, j, k]     = check_large(opa.κ[i, j, k],   large)  ? opa.κ[i, j, k]    : large
                    if has_src
                        opa.src[i, j, k]   = check_large(opa.src[i, j, k], large)  ? opa.src[i, j, k]  : large
                    end
                end
            end
        end
    end
end

set_limits!(aos, opa=nothing; small=1e-30, large=1e30) = begin
    set_small!(aos, opa, small)
    set_large!(aos, opa, large)
end

"""
    smoothAccumulate(eos)

Smooth maxima in EoS lnT-lnEi relation that cause problems when interpolating 
by reverse accumulating the minima and fitting a smoothing spline (optional).
This should leave monotonic arrays unchanged and make others monotonically increasing.
"""
smoothAccumulate!(aos; spline=true) = begin
    @assert !is_internal_energy(aos)

    eos = aos.eos

    x = copy(eos.lnT)
	y = similar(eos.lnT)

    for i in axes(eos.lnEi, 2)
        y .= reverse(eos.lnEi[:, i])
        y .= reverse(accumulate(min, y))

        eos.lnEi[:, i] .= if spline
            ip = Interpolations.extrapolate(
                Interpolations.interpolate(
                    Interpolations.deduplicate_knots!(
                        x, 
                        move_knots=true
                    ), 
                    y, 
                    Interpolations.SteffenMonotonicInterpolation()
                ), 
                Interpolations.Flat()
            )
	        ip.(x) 
        else
            y
        end
    end

    aos
end
smoothAccumulate!(eos::EoSTable, args...; kwargs...) = smoothAccumulate!(@axed(eos), args...; kwargs...)


check(v, small) = !isnan(v) & (v >= small)
check_small = check
check_large(v, large) = !isnan(v) & (v < large)
