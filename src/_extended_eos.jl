# ============================================================================
# Extended EoS -- Flexible extension of the EoS for arbitrary quantities
# ============================================================================
@kwdef struct ExtendedEoS{E<:AxedEoS}
	eos::E
	extensions::Dict = Dict()
end

@kwdef struct ExtendedOpacity
	opa 
	extensions::Dict = Dict()
	binned::Bool = all(diff(opa.λ) .== 1) ? true : false
	weights = binned ? ones(eltype(opa.κ), length(opa.λ)) : ω_midpoint(opa)
end

# ============================================================================
# Access functions
# ============================================================================

opacity(opa::ExtendedOpacity) = opa.opa
wavelength(opa::ExtendedOpacity) = wavelength(opa.opa)
extended(opa::RegularOpacityTable) = ExtendedOpacity(opa=opa)
extended(eos::RegularEoSTable) = ExtendedEoS(eos=@axed(eos))
extended(eos::AxedEoS) = ExtendedEoS(eos=eos)

get_data(eos::ExtendedEoS, var::Symbol) = begin
    Typ = eltype(eos.eos.eos.lnPg)
    if var in fieldnames(typeof(eos.eos.eos))
        return getfield(eos.eos.eos, var)::AbstractArray{Typ}
    elseif var in keys(eos.extensions)
        return eos.extensions[var]::AbstractArray{Typ}
    else
        error("Variable $var not found in EoS tables")
    end
end

get_data(eos::ExtendedEoS, opa::ExtendedOpacity, var::Symbol) = begin
    Typ = eltype(eos.eos.eos.lnPg)
    if var in fieldnames(typeof(eos.eos.eos))
        return getfield(eos.eos.eos, var)::AbstractArray{Typ}
    elseif var in keys(eos.extensions)
        return eos.extensions[var]::AbstractArray{Typ}
    elseif var in fieldnames(typeof(opa.opa))
        return getfield(opa.opa, var)::AbstractArray{Typ}
    elseif var in keys(opa.extensions)
        return opa.extensions[var]::AbstractArray{Typ}
    else
        error("Variable $var not found in tables")
    end
end

# ============================================================================
# Extend general lookup interface
# ============================================================================

function extended_lookup(eos_ext::ExtendedEoS, what, pressure_var, energy_var; maxiter=200, tol=1e-8)
    eos = eos_ext.eos
    if what == :lnRho
        rho_min = minimum(eos.eos.lnRho)
        rho_max = maximum(eos.eos.lnRho)
        rho_mid = Ref(1.0)
        for i in 1:maxiter
            rho_mid[] = (rho_min + rho_max)/2
            P_calc = lookup(eos.eos, :lnPg, rho_mid[], energy_var)
			if abs(rho_max - rho_min) <= tol
	            break
	        end
            if P_calc < pressure_var
                rho_min = rho_mid[]
            else
                rho_max = rho_mid[]
            end
        end
        rho_mid[]
	elseif what in [:lnEi, :lnPg, :lnRoss, :lnEi]
        lookup(eos.eos, what, pressure_var, energy_var)
	elseif what in keys(eos_ext.extensions)
		x = eos.eos.lnT
		y = eos.eos.lnRho
		z = eos_ext.extensions[what]
		xy = TSO.linear_interpolation((x, y), z, extrapolation_bc=Line())
		lookup(TSO.UniformLookup(xy), pressure_var, energy_var)
	else
		error("Given $(what) is not included in the extenden EoS.")
    end
end

function extended_lookup(eos, opa::ExtendedOpacity, what, density_var, energy_var, args...)
    if what in fieldnames(RegularOpacityTable)
        lookup(eos, opa, what, density_var, energy_var)
	elseif what in keys(opa.extensions)
		second_axis = EnergyAxis(eos).unirange
        rho_axis    = DensityAxis(eos).unirange
        ip = linear_interpolation((second_axis, rho_axis), view(opa.extensions[what], :, :, args...), extrapolation_bc=Line())
		u = UniformLookup((x1, x2) -> ip(x1, x2))
		lookup(u, density_var, energy_var)
	else
		error("Given $(what) is not included in the extenden Opacity.")
    end
end

lookup(eos::ExtendedEoS, opa::RegularOpacityTable, args...; kwargs...) = lookup(eos.eos, opa, args...; kwargs...)
lookup(eos::ExtendedEoS, args...; kwargs...) = extended_lookup(eos, args...; kwargs...)
lookup(eos, opa::ExtendedOpacity, args...; kwargs...) = extended_lookup(eos, opa, args...; kwargs...)

# ============================================================================
# Bilinear interpolation coefficients
# ============================================================================

global LOG_INTERPOLATED_VARS = [:κ, :κ_ross, :src]
struct InterpCoefs{T<:AbstractFloat}
    idx::Int        
    w_low::T  
    w_high::T 
end

_is_log_interpolated(var::Symbol) = var in LOG_INTERPOLATED_VARS
add_log_interpolated!(var::Symbol, to=LOG_INTERPOLATED_VARS) = push!(to, var)

weights(eos::ExtendedEoS, lnrho::AbstractArray, lnT::AbstractArray) = begin
    #grid_T = energy_variable(eos.eos) == :lnT ? eos.eos.eos.lnT : eos.eos.eos.lnEi # Assume lnT input represents the energy variable of the table
    grid_T = EnergyAxis(eos.eos).values
    grid_Rho = DensityAxis(eos.eos).values
    
    coefs_T = weights_axis(lnT, grid_T)
    coefs_Rho = weights_axis(lnrho, grid_Rho)
    
    return coefs_Rho, coefs_T
end

weights_axis(x, grid) = begin
    T = eltype(x)
    coefs = Array{InterpCoefs{T}}(undef, size(x))
    
    @inbounds for j in eachindex(x)
        coefs[j] = linear_interpolation_weights(grid, x[j])
    end
    return coefs
end

function linear_interpolation_weights(grid_nodes, val)
    #T = eltype(grid_nodes)
    T = eltype(val)
    val_T = Base.convert(T, val)
    
    if val_T <= first(grid_nodes)
        return InterpCoefs{T}(1, one(T), zero(T))
    elseif val_T >= last(grid_nodes)
        return InterpCoefs{T}(length(grid_nodes)-1, zero(T), one(T))
    end
    
    i = searchsortedlast(grid_nodes, val_T)
    
    x0 = grid_nodes[i]
    x1 = grid_nodes[i+1]
    w_high = (val_T - x0) / (x1 - x0)
    w_low  = one(T) - w_high
    
    return InterpCoefs{T}(i, w_low, w_high)
end
    
# ============================================================================
# Fast lookup function for multiple quantities
# ============================================================================

sample(eos::RegularEoSTable, args...; kwargs...) =  sample(extended(eos), args...; kwargs...)
sample(eos::AxedEoS, args...; kwargs...) = sample(extended(eos), args...; kwargs...)
sample(eos::ExtendedEoS, opa::RegularOpacityTable, args...; kwargs...) = sample(eos, extended(opa), args...; kwargs...)

sample(eos::ExtendedEoS, variables, coord1::Number, coord2::Number, args...) = begin
    res = sample(eos, variables, [coord1], [coord2], args...)
    return map(r -> ndims(r) == 1 ? r[1] : vec(r[:, 1]), res)
end
sample(eos::ExtendedEoS, opa::ExtendedOpacity, variables, coord1::Number, coord2::Number, args...) = begin
    res = sample(eos, opa, variables, [coord1], [coord2], args...)
    return map(r -> ndims(r) == 1 ? r[1] : vec(r[:, 1]), res)
end

function sample(eos::ExtendedEoS, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...) where {T<:AbstractFloat}
    lnrho = coord1
    lnT   = coord2
    
    results = map(var -> similar(lnrho, eltype(lnrho)), variables)
    
    sample!(results, eos, variables, lnrho, lnT, args...)
    results
end

function sample(eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...) where {T<:AbstractFloat}
    lnrho = coord1
    lnT   = coord2

    results = map(variables) do var
        if var == :lnRho || (var == energy_variable(eos.eos))
            similar(lnrho, eltype(lnrho))
        else
            data = get_data(eos, opa, var)
            if ndims(data) == 2
                similar(lnrho, eltype(lnrho))
            else
                if length(args) > 0
                    k = args[1]
                    if k isa Integer
                        similar(lnrho, eltype(lnrho))
                    else
                        similar(lnrho, eltype(lnrho), (length(k), size(lnrho)...))
                    end
                else
                    similar(lnrho, eltype(lnrho), (size(data, 3), size(lnrho)...))
                end
            end
        end
    end
    sample!(results, eos, opa, variables, lnrho, lnT, args...)
    results
end

sample!(results::Union{Tuple, AbstractArray}, eos::RegularEoSTable, args...; kwargs...) = sample!(results, extended(eos), args...; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::AxedEoS, args...; kwargs...) = sample!(results, extended(eos), args...; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::RegularOpacityTable, args...; kwargs...) = sample!(results, eos, extended(opa), args...; kwargs...)

sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...) where {T<:AbstractFloat} = begin
    lnrho, lnT = _invert_axis(eos, variables, coord1, coord2)
    coefs_Rho, coefs_T = weights(eos, lnrho, lnT)
    sample!(results, eos, variables, coefs_Rho, coefs_T, args...)
    
    results
end
sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...) where {T<:AbstractFloat} = begin
    lnrho, lnT = _invert_axis(eos, variables, coord1, coord2)
    coefs_Rho, coefs_T = weights(eos, lnrho, lnT)
    sample!(results, eos, opa, variables, coefs_Rho, coefs_T, args...)
    results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        if var == :lnRho
            grid = DensityAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_Rho, false)
        elseif var == energy_variable(eos.eos)
            grid = EnergyAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_T, false)
        else
            data = get_data(eos, var)
            is_log = _is_log_interpolated(var)
            _interpolate_var!(results[i], data, coefs_Rho, coefs_T, is_log)
        end
    end
    return results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        if var == :lnRho
            grid = DensityAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_Rho, false)
        elseif var == energy_variable(eos.eos)
            grid = EnergyAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_T, false)
        else
            data = get_data(eos, opa, var)
            is_log = _is_log_interpolated(var)
            _interpolate_var!(results[i], data, coefs_Rho, coefs_T, is_log)
        end
    end
    
    return results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}, k::Union{Integer, AbstractUnitRange}) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        if var == :lnRho
            grid = DensityAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_Rho, false)
        elseif var == energy_variable(eos.eos)
            grid = EnergyAxis(eos.eos).values
            _interpolate_var!(results[i], grid, coefs_T, false)
        else
            data = get_data(eos, opa, var)
            is_log = _is_log_interpolated(var)
            _interpolate_var!(results[i], data, coefs_Rho, coefs_T, k, is_log)
        end
    end
    
    return results
end

# ============================================================================
# Inverse lookup for EoS coordinates
# ============================================================================

_invert_axis(eos::ExtendedEoS, variables, coord1::AbstractArray{T}, coord2::AbstractArray{T}) where {T<:AbstractFloat} = begin
    has_target_rho = :lnRho in variables
    grid_energy_var = energy_variable(eos.eos)
    dependent_energy_var = dependent_energy_variable(eos.eos)
    has_target_energy = grid_energy_var in variables
    
    lnrho, lnT = if has_target_rho
        # coord1 is lnPg, coord2 is grid_energy
        _invert(eos, :lnRho, :lnPg, coord1, coord2), coord2
    elseif has_target_energy
        # coord1 is lnRho, coord2 is known_energy (the non-grid one)
        coord1, _invert(eos, grid_energy_var, dependent_energy_var, coord2, coord1)
    else
        coord1, coord2
    end
    return lnrho, lnT
end

function _invert(eos::ExtendedEoS, target_var::Symbol, known_var::Symbol, target_array::AbstractArray{T}, known_array::AbstractArray{T}; maxiter=200, tol=1e-6) where {T}
    result_array = similar(target_array)
    grid = getfield(eos.eos.eos, target_var)
    val_min_global = minimum(grid)
    val_max_global = maximum(grid)
    
    calc_arr = zeros(T, 1)
    results = (calc_arr,)
    mid_arr = zeros(T, 1)
    known_val_arr = zeros(T, 1)
    
    is_finding_rho = target_var == :lnRho
    
    @inline function eval_coord(coord_guess, coefs_Rho, coefs_T)
        mid_arr[1] = coord_guess
        
        if is_finding_rho
            coefs_Rho .= weights_axis(mid_arr, grid)
        else
            coefs_T .= weights_axis(mid_arr, grid)
        end
        
        sample!(results, eos, (known_var,), coefs_Rho, coefs_T)
        return calc_arr[1]
    end

    is_invalid(f_min, f_max) = (sign(f_min) == sign(f_max)) || (is_finding_rho && (f_max < f_min))

    calc_fallback(v_min, v_max, f_min, f_max, target_val) = begin
        if is_finding_rho
            lnP_val = target_val
            if energy_variable(eos.eos) == :lnT
                lnT_val = known_val_arr[1]
            else
                a_rad = 7.5657e-15 # erg cm^-3 K^-4
                lnT_val = 0.25 * (lnP_val + log(3.0) - log(a_rad))
            end
            
            mu_guess = 0.607
            return lnP_val + log(mu_guess * m_u) - log(KBoltzmann) - lnT_val
        else
            return abs(f_min) < abs(f_max) ? v_min : v_max
        end
    end
    
    @inbounds for j in eachindex(target_array)
        target_val = target_array[j]
        known_val_arr[1] = known_array[j]
        v_mid = (val_min_global + val_max_global) / 2
        mid_arr[1] = v_mid
        
        coefs_Rho, coefs_T = if is_finding_rho
            weights(eos, mid_arr, known_val_arr)
        else
            weights(eos, known_val_arr, mid_arr)
        end
        
        v_ideal = calc_fallback(val_min_global, val_max_global, 0.0, 0.0, target_val)
        
        v_min = val_min_global
        v_max = val_max_global
        
        if is_finding_rho
            v_min_tight = max(val_min_global, v_ideal - 2.0)
            v_max_tight = min(val_max_global, v_ideal + 2.0)
            
            f_min_tight = eval_coord(v_min_tight, coefs_Rho, coefs_T) - target_val
            f_max_tight = eval_coord(v_max_tight, coefs_Rho, coefs_T) - target_val
            
            if sign(f_min_tight) != sign(f_max_tight) && f_max_tight > f_min_tight
                v_min = v_min_tight
                v_max = v_max_tight
                f_min = f_min_tight
                f_max = f_max_tight
            else
                f_min = eval_coord(v_min, coefs_Rho, coefs_T) - target_val
                f_max = eval_coord(v_max, coefs_Rho, coefs_T) - target_val
            end
        else
            f_min = eval_coord(v_min, coefs_Rho, coefs_T) - target_val
            f_max = eval_coord(v_max, coefs_Rho, coefs_T) - target_val
        end

        for i in 1:maxiter
            v_mid = (v_min + v_max)/2
            f_mid = eval_coord(v_mid, coefs_Rho, coefs_T) - target_val
            
            if abs(v_max - v_min) <= tol || f_mid == 0
                result_array[j] = v_mid
                break
            end
            
            if is_invalid(f_min, f_max)
                if v_ideal > v_mid
                    v_min = v_min + abs(v_min - v_mid) / 2.0
                    v_max = min(val_max_global, v_max + abs(v_min - v_mid) / 2.0)
                else
                    v_max = v_max - abs(v_max - v_mid) / 2.0
                    v_min = max(val_min_global, v_min - abs(v_max - v_mid) / 2.0)
                end
                
                f_min = eval_coord(v_min, coefs_Rho, coefs_T) - target_val
                f_max = eval_coord(v_max, coefs_Rho, coefs_T) - target_val
            else
                # Regular bisection
                if sign(f_mid) == sign(f_min)
                    v_min = v_mid
                    f_min = f_mid
                else
                    v_max = v_mid
                    f_max = f_mid
                end
            end
            
            if i == maxiter
                result_array[j] = v_mid
            end
        end
    end
    result_array
end

# ============================================================================
# Core interpolation
# ============================================================================

_interpolate_var!(result::AbstractArray{T3}, data::AbstractArray{T}, coefs::AbstractArray{InterpCoefs{T2}}, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat} = begin
    @inbounds for j in 1:length(coefs)
        result[j] = _interpolate_point_2d(data, coefs[j], is_log)
    end
end

_interpolate_var!(result::AbstractArray{T3}, data::AbstractArray{T}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat} = begin
    if ndims(data) == 2
        @inbounds for j in 1:length(coefs_T)
            result[j] = _interpolate_point_2d(data, coefs_Rho[j], coefs_T[j], is_log)
        end
    elseif ndims(data) == 3
        n3 = size(data, 3)
        @inbounds for j in 1:length(coefs_T)
            cr = coefs_Rho[j]
            ct = coefs_T[j]
            
            it, ir = ct.idx, cr.idx
            w00 = ct.w_low  * cr.w_low
            w10 = ct.w_high * cr.w_low
            w01 = ct.w_low  * cr.w_high
            w11 = ct.w_high * cr.w_high
            
            for k in 1:n3
                result[k + (j-1)*n3] = _interpolate_point_3d_fused(data, it, ir, w00, w10, w01, w11, k, is_log)
            end
        end
    end
end

_interpolate_var!(result::AbstractArray{T3}, data::AbstractArray{T}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}, k::Union{Integer, AbstractUnitRange}, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat} = begin
    if ndims(data) == 2
        @inbounds for j in 1:length(coefs_T)
            result[j] = _interpolate_point_2d(data, coefs_Rho[j], coefs_T[j], is_log)
        end
    elseif ndims(data) == 3
        if k isa Integer
            @inbounds for j in 1:length(coefs_T)
                cr = coefs_Rho[j]
                ct = coefs_T[j]
                it, ir = ct.idx, cr.idx
                w00 = ct.w_low  * cr.w_low
                w10 = ct.w_high * cr.w_low
                w01 = ct.w_low  * cr.w_high
                w11 = ct.w_high * cr.w_high
                
                result[j] = _interpolate_point_3d_fused(data, it, ir, w00, w10, w01, w11, k, is_log)
            end
        else
            nk = length(k)
            @inbounds for j in 1:length(coefs_T)
                cr = coefs_Rho[j]
                ct = coefs_T[j]
                it, ir = ct.idx, cr.idx
                w00 = ct.w_low  * cr.w_low
                w10 = ct.w_high * cr.w_low
                w01 = ct.w_low  * cr.w_high
                w11 = ct.w_high * cr.w_high
                
                for (ik, kidx) in enumerate(k)
                    result[ik + (j-1)*nk] = _interpolate_point_3d_fused(data, it, ir, w00, w10, w01, w11, kidx, is_log)
                end
            end
        end
    end
end

@inline _interpolate_point_2d(data::AbstractArray{T}, coefs::InterpCoefs{T2}, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat} = begin
    it = coefs.idx
    w00 = coefs.w_low
    w10 = coefs.w_high
    
    if is_log
        val = w00 * log(data[it]) +
              w10 * log(data[it+1])
        return exp(val)
    else
        return w00 * data[it] +
               w10 * data[it+1]
    end
end

@inline function _interpolate_point_2d(data::AbstractArray{T}, cr::InterpCoefs{T2}, ct::InterpCoefs{T2}, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat}
    it, ir = ct.idx, cr.idx
    w00 = ct.w_low  * cr.w_low
    w10 = ct.w_high * cr.w_low
    w01 = ct.w_low  * cr.w_high
    w11 = ct.w_high * cr.w_high
    
    if is_log
        val = w00 * log(data[it,   ir  ]) +
              w10 * log(data[it+1, ir  ]) +
              w01 * log(data[it,   ir+1]) +
              w11 * log(data[it+1, ir+1])
        return exp(val)
    else
        return w00 * data[it,   ir  ] +
               w10 * data[it+1, ir  ] +
               w01 * data[it,   ir+1] +
               w11 * data[it+1, ir+1]
    end
end

@inline function _interpolate_point_3d_fused(data::AbstractArray{T}, it::Int, ir::Int, w00::T2, w10::T2, w01::T2, w11::T2, k::Int, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat}
    if is_log
        val = w00 * log(data[it,   ir,   k]) +
              w10 * log(data[it+1, ir,   k]) +
              w01 * log(data[it,   ir+1, k]) +
              w11 * log(data[it+1, ir+1, k])
        return exp(val)
    else
        return w00 * data[it,   ir,   k] +
               w10 * data[it+1, ir,   k] +
               w01 * data[it,   ir+1, k] +
               w11 * data[it+1, ir+1, k]
    end
end

# ============================================================================
# Documentation
# ============================================================================

"""
    sample!(results, eos, opa, variables, lnrho, lnT, [k])
    sample!(results, eos, opa, variables, coefs_Rho, coefs_T, [k])

In-place, mutating version of `sample`.

This function is optimized for repeatedly evaluating tables without memory allocations. 
By manually passing the precomputed interpolation weights, you avoid binary tree search grid lookups entirely!

# Examples
```julia
# 1) Pre-allocate output arrays matching your array variable
results = (similar(lnrho), similar(lnrho))

# 2) Pre-compute interpolation coefficients for a fixed grid of coordinates
w_Rho, w_T = weights(eos, lnrho, lnT)

# 3) Repeatedly evaluate over different wavelength slices without any allocations or indexing algorithms!
for k in 1:100
    sample!(results, eos, opa, (:κ, :src), w_Rho, w_T, k)
    # ... use results[1] and results[2] matrix blocks
end
```
"""
sample!

"""
    sample(eos, opa, variables, lnrho, lnT, [k])

Perform highly optimized interpolation on extended EOS and opacity tables.

# Arguments
- `eos`: An `ExtendedEoS` table.
- `opa`: An `ExtendedOpacity` table.
- `variables`: A collection (e.g., `Tuple` or `AbstractArray`) of symbols measuring what to interpolate (e.g., `[:κ, :κ_ross]`).
- `lnrho`, `lnT`: Logarithm values of the density and temperature. Can be a scalar `Number` or an `AbstractArray{<:AbstractFloat}`.
- `k` (optional): An `Integer` index or `AbstractUnitRange` of indices to slice the 3rd dimension (e.g., wavelength) of 3D opacity arrays. If omitted, interpolates all wavelengths.

# Examples
```julia
# Scalar point lookup
κ, κ_ross = sample(eos, opa, (:κ, :κ_ross), -10.0, 8.5)

# Array lookup
κ_arr, = sample(eos, opa, [:κ], lnrho_array, lnt_array)

# Array lookup sliced over a wavelength chunk
κ_slice, = sample(eos, opa, [:κ], lnrho_array, lnt_array, 1:50)
```
"""
sample

# ============================================================================
# Ideal gas entropy
# ============================================================================

function ideal_entropy(T, rho, mu)    
    m_particle = mu * m_u
    R_spec = KBoltzmann / m_particle
    
    factor = (2.0 * pi * m_particle * KBoltzmann * T) / (HPlanck^2)
    n_quantum = factor^1.5
    n = rho / m_particle
    
    # Sackur-Tetrode Equation
    R_spec * (log(n_quantum / n) + 2.5)
end

# ============================================================================
# Table gradients for pre-computation
# ============================================================================

function gradients!(eos, opa_extended::ExtendedOpacity)
	aos = @axed eos
	@assert !is_internal_energy(aos)
	opa = opa_extended.opa

	nT, nRho = size(eos)
	dS_dT = similar(opa.src)

	Threads.@threads for k in axes(opa.src, 3)
		@inbounds for j in eachindex(eos.lnRho)
			@inbounds for i in eachindex(eos.lnT)
				dS, dT = if (i>1) &&(i<nT)
					log(opa.src[i+1, j, k]) - log(opa.src[i-1, j, k]), 
					eos.lnT[i+1] - eos.lnT[i-1]
				elseif i==1
					log(opa.src[i+1, j, k]) - log(opa.src[i, j, k]), 
					eos.lnT[i+1] - eos.lnT[i]
				elseif i==nT
					log(opa.src[i, j, k]) - log(opa.src[i-1, j, k]), 
					eos.lnT[i] - eos.lnT[i-1]
				end
				dS_dT[i, j, k] = exp(log(opa.src[i, j, k]) - eos.lnT[i]) * dS / dT
			end
		end
	end

	opa_extended.extensions[:dS_dT] = dS_dT

	dS_dT
end

"""
	add_gradients!(eos)

Compute χₜ, χᵨ, and cᵥ for the given EoS. 
"""
function gradients!(eos_extended::ExtendedEoS)
	eos = eos_extended.eos.eos
	aos = eos_extended.eos
	@assert !is_internal_energy(aos)

	nT, nRho = size(eos)
	χₜ = similar(eos.lnPg)
	χᵨ = similar(eos.lnPg)
	cᵥ = similar(eos.lnPg)
	dlnRoss_dlnT = similar(eos.lnPg)
	dlnRoss_dlnRho = similar(eos.lnPg)

	@inbounds for j in eachindex(eos.lnRho)
		@inbounds for i in eachindex(eos.lnT)
			dE, dT, dP, dRoss = if (i>1) &&(i<nT)
				eos.lnEi[i+1, j] - eos.lnEi[i-1, j], 
				eos.lnT[i+1    ] - eos.lnT[i-1],
				eos.lnPg[i+1, j] - eos.lnPg[i-1, j],
				eos.lnRoss[i+1, j] - eos.lnRoss[i-1, j]
			elseif i==1
				eos.lnEi[i+1, j] - eos.lnEi[i, j], 
				eos.lnT[i+1]     - eos.lnT[i],
				eos.lnPg[i+1, j] - eos.lnPg[i, j],
				eos.lnRoss[i+1, j] - eos.lnRoss[i, j]
			elseif i==nT
				eos.lnEi[i, j] - eos.lnEi[i-1, j], 
				eos.lnT[i]     - eos.lnT[i-1],
				eos.lnPg[i, j] - eos.lnPg[i-1, j],
				eos.lnRoss[i, j] - eos.lnRoss[i-1, j]
			end
			cᵥ[i, j] = max(dE/dT * exp(eos.lnEi[i, j] - eos.lnT[i]), 1e-12)
			χₜ[i, j] = max(dP/dT, 1e-12)
			dlnRoss_dlnT[i, j] = dRoss / dT
	
			dP, dR, dRoss = if (j>1) &&(j<nRho)
				eos.lnRho[j+1] - eos.lnRho[j-1],
				eos.lnPg[i, j+1] - eos.lnPg[i, j-1],
				eos.lnRoss[i, j+1] - eos.lnRoss[i, j-1]
			elseif j==1
				eos.lnRho[j+1] - eos.lnRho[j],
				eos.lnPg[i, j+1] - eos.lnPg[i, j],
				eos.lnRoss[i, j+1] - eos.lnRoss[i, j]
			elseif j==nRho
				eos.lnRho[j] - eos.lnRho[j-1],
				eos.lnPg[i, j] - eos.lnPg[i, j-1],
				eos.lnRoss[i, j] - eos.lnRoss[i, j-1]
			end
			χᵨ[i, j] = max(dP/dR, 1e-12)
			dlnRoss_dlnRho[i, j] = dRoss / dR
			dlnRoss_dlnT[i, j] = dlnRoss_dlnT[i, j] - χₜ[i, j] / χᵨ[i, j]*dlnRoss_dlnRho[i, j]
		end
	end

	eos_extended.extensions[:χₜ] = χₜ
	eos_extended.extensions[:χᵨ] = χᵨ
	eos_extended.extensions[:cᵥ] = cᵥ
	eos_extended.extensions[:dlnκ_dlnT] = dlnRoss_dlnT

	χₜ, χᵨ, cᵥ
end

# ============================================================================
# Thermodynamic quantities for pre-computation
# ============================================================================

"""
    add_thermodynamics!(eos_extended)

Add gradients and thermodynamic quantities to the EoS and allow for interpolation.
"""
function add_thermodynamics!(eos_extended::ExtendedEoS)
	eos = eos_extended.eos.eos
	χₜ, χᵨ, cᵥ = gradients!(eos_extended)

	# expansion coefficient dlnrho / dlnT |_P = χₜ/χᵨ
	Q = χₜ ./ χᵨ

	# specific heat cp = cv + P/(rho T) χₜ²/χᵨ
	cₚ = similar(cᵥ)
	∇ₐ = similar(cᵥ)
	μ = similar(cᵥ)
	S = zeros(size(cᵥ))
	for j in eachindex(eos.lnRho)
		for i in eachindex(eos.lnT)
			P = exp(eos.lnPg[i, j])
			ρ = exp(eos.lnRho[j])
			T = exp(eos.lnT[i])
			cₚ[i, j] = min(1e12, max(1e-12, cᵥ[i, j] + P/(ρ*T) * χₜ[i, j]^2 / χᵨ[i, j]))
			∇ₐ[i, j] = min(1, max(1e-12, P/(ρ*T) .* Q[i, j] ./ cₚ[i, j]))
			μ[i, j] = min(1e12, max(1e-12, ρ * KBoltzmann * T / (P * m_u)))
		end
	end
	eos_extended.extensions[:Q] = Q
	eos_extended.extensions[:cₚ] = cₚ
	eos_extended.extensions[:∇ₐ] = ∇ₐ
	eos_extended.extensions[:μ] = μ

	# compute the entropy by starting with the ideal gas value and then integrating
	# along isotherms. First, we integrate along isobar (constant rho)
	S[1, 1] = ideal_entropy(exp(eos.lnT[1]), exp(eos.lnRho[1]), μ[1, 1])
	for i in 2:size(S, 1)
		dE = eos.lnEi[i, 1] - eos.lnEi[i-1, 1]
		T = (exp(eos.lnT[i]) +  exp(eos.lnT[i-1])) / 2
		E = (exp(eos.lnEi[i, 1]) +  exp(eos.lnEi[i-1, 1])) / 2
		
		dS = E / T * dE	
		S[i, 1] = S[i-1, 1] + dS
	end

	# and then along constant T
	T = exp.(eos.lnT)
	for j in 2:size(S, 2)
		dE = eos.lnEi[:, j] - eos.lnEi[:, j-1]
		dRho = eos.lnRho[j] - eos.lnRho[j-1]
								   
		rho = (exp(eos.lnRho[j]) + exp(eos.lnRho[j-1])) / 2
		P = (exp.(eos.lnPg[:, j]) + exp.(eos.lnPg[:, j-1])) / 2
		E = (exp.(eos.lnEi[:, j]) + exp.(eos.lnEi[:, j-1])) / 2

		t1 = E ./ T .* dE
		t2 = P ./ (rho .* T) .* dRho
		
		dS = t1 - t2
		S[:, j] = S[:, j-1] .+ dS
	end
	eos_extended.extensions[:S] = S

	nothing
end