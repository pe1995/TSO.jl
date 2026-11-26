# ============================================================================
# Sampling of EquationOfState tables
# ============================================================================

# ============================================================================
# Bilinear interpolation coefficients
# ============================================================================

global LOG_INTERPOLATED_VARS = [:κ, :src]
struct InterpCoefs{T<:AbstractFloat}
    idx::Int
    w_low::T
    w_high::T
end

_is_log_interpolated(var::Symbol) = var in LOG_INTERPOLATED_VARS
add_log_interpolated!(var::Symbol, to=LOG_INTERPOLATED_VARS) = push!(to, var)

weights(eos::EquationOfState, lnrho::AbstractArray, lnT::AbstractArray; extrapolate=false) = begin
    grid_T = eos.energy_axis.values
    grid_Rho = eos.density_axis.values

    coefs_T = weights_axis(lnT, grid_T, extrapolate)
    coefs_Rho = weights_axis(lnrho, grid_Rho, extrapolate)

    return coefs_Rho, coefs_T
end

"""
    weights(eos; coordinates...)

Keyword version, accepting the same named (and possibly linear) coordinates as [`sample`](@ref),
e.g. `weights(eos, rho=rho, T=T)`.
"""
weights(eos::EquationOfState; extrapolate=false, kwargs...) = begin
    given, _ = _coord_kwargs(kwargs)
    weights(eos, _grid_coords(eos, given; extrapolate=extrapolate)...; extrapolate=extrapolate)
end

weights_axis(x, grid, extrapolate=false) = begin
    T = eltype(x)
    coefs = Array{InterpCoefs{T}}(undef, size(x))
    weights_axis!(coefs, x, grid, extrapolate)
end

weights_axis!(coefs, x, grid, extrapolate=false) = begin
    @inbounds for j in eachindex(x)
        coefs[j] = linear_interpolation_weights(grid, x[j], extrapolate)
    end
    return coefs
end

"""
    weights!(coefs_Rho, coefs_T, eos, lnrho, lnT)

In-place [`weights`](@ref): fill caller-owned coefficient arrays instead of
allocating new ones.

For callers that evaluate the tables repeatedly on arrays of the same size -- a
hydrodynamics solver sampling the equation of state once per cell per timestep,
say -- where the two coefficient arrays would otherwise be allocated and
discarded on every call.
"""
weights!(coefs_Rho::AbstractArray{InterpCoefs{T}}, coefs_T::AbstractArray{InterpCoefs{T}},
         eos::EquationOfState, lnrho::AbstractArray, lnT::AbstractArray; extrapolate=false) where {T<:AbstractFloat} = begin
    weights_axis!(coefs_T, lnT, eos.energy_axis.values, extrapolate)
    weights_axis!(coefs_Rho, lnrho, eos.density_axis.values, extrapolate)
    return coefs_Rho, coefs_T
end

function linear_interpolation_weights(grid_nodes, val, extrapolate=false)
    T = eltype(val)
    val_T = Base.convert(T, val)

    if isnan(val_T)
        return InterpCoefs{T}(1, T(NaN), T(NaN))
    end

    n = length(grid_nodes)
    i = if val_T <= first(grid_nodes)
        extrapolate || return InterpCoefs{T}(1, one(T), zero(T))
        1
    elseif val_T >= last(grid_nodes)
        extrapolate || return InterpCoefs{T}(n-1, zero(T), one(T))
        n-1
    else
        clamp(searchsortedlast(grid_nodes, val_T), 1, n-1)
    end

    x0 = grid_nodes[i]
    x1 = grid_nodes[i+1]
    w_high = (val_T - x0) / (x1 - x0)
    w_low  = one(T) - w_high

    return InterpCoefs{T}(i, w_low, w_high)
end

# ============================================================================
# Fast lookup function for multiple quantities
# ============================================================================

get_data(eos::EquationOfState, var::Symbol) = eos[var]

sample(eos::EquationOfState, opa::EquationOfState, args...; kwargs...) = sample(merge(eos, opa), args...; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::EquationOfState, opa::EquationOfState, args...; kwargs...) = sample!(results, merge(eos, opa), args...; kwargs...)

sample(eos::EquationOfState, variables::Union{Tuple, AbstractArray}, coord1::Number, coord2::Number, args...; kwargs...) = begin
    res = sample(eos, variables, [coord1], [coord2], args...; kwargs...)
    map(r -> ndims(r) == 1 ? r[1] : vec(r[:, 1]), res)
end

function sample(eos::EquationOfState, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...; extrapolate=false) where {T<:AbstractFloat}
    results = _allocate_results(eos, variables, coord1, args...)
    sample!(results, eos, variables, coord1, coord2, args...; extrapolate=extrapolate)
    results
end

function _allocate_results(eos::EquationOfState, variables::Union{Tuple, AbstractArray}, lnrho, args...)
    map(variables) do var
        if var == :lnRho || var == energy_variable(eos) || ndims(get_data(eos, var)) == 2
            similar(lnrho, eltype(lnrho))
        elseif length(args) > 0
            k = args[1]
            k isa Integer ? similar(lnrho, eltype(lnrho)) : similar(lnrho, eltype(lnrho), (length(k), size(lnrho)...))
        else
            similar(lnrho, eltype(lnrho), (size(get_data(eos, var), 3), size(lnrho)...))
        end
    end
end

sample!(results::Union{Tuple, AbstractArray}, eos::EquationOfState, variables::Union{Tuple, AbstractArray}, coord1::AbstractArray{T}, coord2::AbstractArray{T}, args...; extrapolate=false) where {T<:AbstractFloat} = begin
    lnrho, lnT = _invert_axis(eos, variables, coord1, coord2; extrapolate=extrapolate)
    coefs_Rho, coefs_T = weights(eos, lnrho, lnT; extrapolate=extrapolate)
    sample!(results, eos, variables, coefs_Rho, coefs_T, args...)
    results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::EquationOfState, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}, k...) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        if var == :lnRho
            _interpolate_var!(results[i], eos.density_axis.values, coefs_Rho, false)
        elseif var == energy_variable(eos)
            _interpolate_var!(results[i], eos.energy_axis.values, coefs_T, false)
        else
            _interpolate_var!(results[i], get_data(eos, var), coefs_Rho, coefs_T, k..., _is_log_interpolated(var))
        end
    end
    results
end

# ============================================================================
# Keyword interface -- name the coordinates instead of ordering them
# ============================================================================

const COORDINATE_KWARGS = Dict(
    :lnRho => :lnRho, :rho => :lnRho, :Rho => :lnRho, :ρ => :lnRho,
    :lnT   => :lnT,   :T   => :lnT,
    :lnEi  => :lnEi,  :Ei  => :lnEi,
    :lnPg  => :lnPg,  :Pg  => :lnPg
)

_needs_log(name::Symbol) = !startswith(String(name), "ln")

sample(eos::EquationOfState, variables::Symbol...; kwargs...) = begin
    res = _sample_kwargs(eos, variables; kwargs...)
    length(variables) == 1 ? first(res) : res
end
sample(eos::EquationOfState, variables::Union{Tuple, AbstractArray}; kwargs...) = _sample_kwargs(eos, variables; kwargs...)

sample!(results::Union{Tuple, AbstractArray}, eos::EquationOfState, variables::Symbol...; kwargs...) = _sample_kwargs!(results, eos, variables; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::EquationOfState, variables::Union{Tuple, AbstractArray}; kwargs...) = _sample_kwargs!(results, eos, variables; kwargs...)

function _sample_kwargs(eos::EquationOfState, variables; k=nothing, coefs_Rho=nothing, coefs_T=nothing, extrapolate=false, kwargs...)
    _check_wavelength_index(eos, k)
    cRho, cT, scalar = _resolve_coefs(eos, coefs_Rho, coefs_T; extrapolate=extrapolate, kwargs...)
    kargs = isnothing(k) ? () : (k,)

    template = similar(cRho, _coef_eltype(cRho))
    results = _allocate_results(eos, variables, template, kargs...)
    sample!(results, eos, variables, cRho, cT, kargs...)

    scalar ? map(r -> ndims(r) == 1 ? r[1] : vec(r[:, 1]), results) : results
end

function _sample_kwargs!(results, eos::EquationOfState, variables; k=nothing, coefs_Rho=nothing, coefs_T=nothing, extrapolate=false, kwargs...)
    _check_wavelength_index(eos, k)
    cRho, cT, _ = _resolve_coefs(eos, coefs_Rho, coefs_T; extrapolate=extrapolate, kwargs...)
    sample!(results, eos, variables, cRho, cT, (isnothing(k) ? () : (k,))...)
    results
end

_check_wavelength_index(eos, k) = (!has_opacity(eos) && !isnothing(k)) && error("The wavelength index `k` needs an opacity table, this table has no opacity.")
_coef_eltype(coefs::AbstractArray{InterpCoefs{T}}) where {T<:AbstractFloat} = T

function _resolve_coefs(eos::EquationOfState, coefs_Rho, coefs_T; extrapolate=false, kwargs...)
    if !isnothing(coefs_Rho) || !isnothing(coefs_T)
        (isnothing(coefs_Rho) || isnothing(coefs_T)) && error("Precomputed weights need both `coefs_Rho` and `coefs_T`.")
        isempty(kwargs) || error("Got coordinates $(Tuple(keys(kwargs))) together with precomputed weights, give either one or the other.")
        return coefs_Rho, coefs_T, false
    end

    given, scalar = _coord_kwargs(kwargs)
    lnrho, lnT = _grid_coords(eos, given; extrapolate=extrapolate)
    cRho, cT = weights(eos, lnrho, lnT; extrapolate=extrapolate)

    cRho, cT, scalar
end

function _coord_kwargs(kwargs)
    given = Dict{Symbol, Any}()
    for (name, value) in kwargs
        haskey(COORDINATE_KWARGS, name) || error("Unknown coordinate keyword `$name`. Available: $(join(sort(String.(collect(keys(COORDINATE_KWARGS)))), ", ")).")
        var = COORDINATE_KWARGS[name]
        haskey(given, var) && error("Coordinate $var was given more than once (via `$name`).")
        given[var] = _needs_log(name) ? log.(value) : value
    end
    isempty(given) && error("No coordinates given. Available: $(join(sort(String.(collect(keys(COORDINATE_KWARGS)))), ", ")).")

    scalar = all(v isa Number for v in values(given))
    T = float(mapreduce(eltype, promote_type, values(given)))
    shape = scalar ? (1,) : size(first(v for v in values(given) if !(v isa Number)))

    for (var, value) in given
        given[var] = if value isa Number
            Base.fill(Base.convert(T, value), shape)
        else
            size(value) == shape || error("Coordinate $var has size $(size(value)), expected $shape.")
            Base.convert(Array{T}, value)
        end
    end

    given, scalar
end

function _grid_coords(eos::EquationOfState, given; extrapolate=false)
    grid_energy_var = energy_variable(eos)
    dependent_energy_var = dependent_energy_variable(eos)

    lnrho = get(given, :lnRho, nothing)
    lnT = get(given, grid_energy_var, nothing)

    if isnothing(lnrho) && !isnothing(lnT) && haskey(given, :lnPg)
        lnrho = _invert(eos, :lnRho, :lnPg, given[:lnPg], lnT; extrapolate=extrapolate)
    end
    if isnothing(lnT) && !isnothing(lnrho) && haskey(given, dependent_energy_var)
        lnT = _invert(eos, grid_energy_var, dependent_energy_var, given[dependent_energy_var], lnrho; extrapolate=extrapolate)
    end

    if isnothing(lnrho) || isnothing(lnT)
        error("Can not construct (:lnRho, $(grid_energy_var)) from the given coordinates $(Tuple(keys(given))). "*
              "Give one of (lnRho, rho, ρ) or (lnPg, Pg), together with one of ($(grid_energy_var), $(dependent_energy_var)) or their linear counterparts.")
    end

    lnrho, lnT
end

# ============================================================================
# Inverse lookup for EoS coordinates
# ============================================================================

_invert_axis(eos::EquationOfState, variables, coord1::AbstractArray{T}, coord2::AbstractArray{T}; extrapolate=false) where {T<:AbstractFloat} = begin
    has_target_rho = :lnRho in variables
    grid_energy_var = energy_variable(eos)
    dependent_energy_var = dependent_energy_variable(eos)
    has_target_energy = grid_energy_var in variables

    lnrho, lnT = if has_target_rho
        _invert(eos, :lnRho, :lnPg, coord1, coord2; extrapolate=extrapolate), coord2
    elseif has_target_energy
        coord1, _invert(eos, grid_energy_var, dependent_energy_var, coord2, coord1; extrapolate=extrapolate)
    else
        coord1, coord2
    end
    return lnrho, lnT
end

function _invert(eos::EquationOfState, target_var::Symbol, known_var::Symbol, target_array::AbstractArray{T}, known_array::AbstractArray{T}; maxiter=200, tol=1e-6, extrapolate=false) where {T}
    result_array = similar(target_array)
    grid = eos[target_var]
    pad = extrapolate ? 0.2 * (maximum(grid) - minimum(grid)) : zero(eltype(grid))
    val_min_global = minimum(grid) - pad
    val_max_global = maximum(grid) + pad

    calc_arr = zeros(T, 1)
    results = (calc_arr,)
    mid_arr = zeros(T, 1)
    known_val_arr = zeros(T, 1)

    is_finding_rho = target_var == :lnRho

    @inline function eval_coord(coord_guess, coefs_Rho, coefs_T)
        mid_arr[1] = coord_guess

        if is_finding_rho
            coefs_Rho .= weights_axis(mid_arr, grid, extrapolate)
        else
            coefs_T .= weights_axis(mid_arr, grid, extrapolate)
        end

        sample!(results, eos, (known_var,), coefs_Rho, coefs_T)
        return calc_arr[1]
    end

    is_invalid(f_min, f_max) = (sign(f_min) == sign(f_max)) || (is_finding_rho && (f_max < f_min))

    calc_fallback(v_min, v_max, f_min, f_max, target_val) = begin
        if is_finding_rho
            lnP_val = target_val
            if energy_variable(eos) == :lnT
                lnT_val = known_val_arr[1]
            else
                a_rad = 7.5657e-15
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
            weights(eos, mid_arr, known_val_arr; extrapolate=extrapolate)
        else
            weights(eos, known_val_arr, mid_arr; extrapolate=extrapolate)
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
    sample!(results, eos, [opa], variables...; coordinates..., [k=k])
    sample!(results, eos, [opa], variables...; coefs_Rho=w_Rho, coefs_T=w_T, [k=k])

In-place, mutating version of `sample`.

This function is optimized for repeatedly evaluating tables without memory allocations.
By manually passing the precomputed interpolation weights, you avoid binary tree search grid lookups entirely!

The keyword form accepts the same named coordinates as [`sample`](@ref), or -- for the
allocation free loop -- the precomputed weights as `coefs_Rho` + `coefs_T`. Coordinates and
weights are mutually exclusive.

# Examples
```julia
# 1) Pre-allocate output arrays matching your array variable
results = (similar(lnrho), similar(lnrho))

# 2) Pre-compute interpolation coefficients for a fixed grid of coordinates
w_Rho, w_T = weights(eos, lnrho, lnT)   # or: weights(eos, rho=rho, T=T)

# 3) Repeatedly evaluate over different wavelength slices without any allocations or indexing algorithms!
for k in 1:100
    sample!(results, eos, opa, (:κ, :src), w_Rho, w_T, k)
    # ... use results[1] and results[2] matrix blocks
end

# The same loop with named arguments
for k in 1:100
    sample!(results, eos, opa, :κ, :src, coefs_Rho=w_Rho, coefs_T=w_T, k=k)
end
```
"""
sample!

"""
    sample(eos, opa, variables, lnrho, lnT, [k])
    sample(eos, [opa], variables...; coordinates..., [k=k])

Perform highly optimized interpolation on extended EOS and opacity tables.

# Arguments
- `eos`: An `EquationOfState`, with its opacity if opacities are sampled.
- `opa`: (optional) A second `EquationOfState` holding the opacity, which is merged with `eos`.
- `variables`: A collection (e.g., `Tuple` or `AbstractArray`) of symbols measuring what to interpolate (e.g., `[:κ, :src]`). In the keyword form the variables may also be passed as plain arguments (e.g. `:κ, :src`).
- `lnrho`, `lnT`: Logarithm values of the density and temperature. Can be a scalar `Number` or an `AbstractArray{<:AbstractFloat}`.
- `k` (optional): An `Integer` index or `AbstractUnitRange` of indices to slice the 3rd dimension (e.g., wavelength) of 3D opacity arrays. If omitted, interpolates all wavelengths.

# Coordinate keywords
Instead of relying on the (variable dependent) order of `coord1, coord2`, the coordinates may
be named: `lnRho`, `lnT`, `lnEi`, `lnPg`, or their linear counterparts `rho` (`ρ`), `T`, `Ei`, `Pg`,
which are converted to their logarithm automatically. Any pair that pins down the two grid axes
of the table works, whatever is missing is inverted from the table, e.g. giving `lnPg` + `lnT`
inverts for the density. Because the coordinates are named, a variable that is both requested
and given (e.g. `sample(eos, :lnPg, :lnT, lnRho=..., lnT=...)`) is simply passed through.
Scalars and arrays may be mixed, a scalar is then broadcast to the shape of the array.
When exactly one variable is requested as a plain argument, its result is returned unwrapped.

# Examples
```julia
# Scalar point lookup
κ, S = sample(eos, opa, (:κ, :src), -10.0, 8.5)

# Array lookup
κ_arr, = sample(eos, opa, [:κ], lnrho_array, lnt_array)

# Array lookup sliced over a wavelength chunk
κ_slice, = sample(eos, opa, [:κ], lnrho_array, lnt_array, 1:50)

# The same, with named coordinates
lnPg = sample(eos, :lnPg, lnRho=lnrho_array, lnT=lnt_array)
lnPg, lnNe = sample(eos, :lnPg, :lnNe, rho=rho_array, T=T_array)
κ_slice = sample(eos, opa, :κ, lnRho=lnrho_array, lnT=lnt_array, k=1:50)
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

    R_spec * (log(n_quantum / n) + 2.5)
end

# ============================================================================
# Table gradients for pre-computation
# ============================================================================

"""
    source_gradient!(tab)

Temperature derivative of the source function, stored as `:dS_dT` in the opacity of `tab`.
"""
function source_gradient!(tab::EquationOfState)
	@assert !is_internal_energy(tab)
	src = tab[:src]
	lnT = tab[:lnT]

	nT, nRho = size(tab)
	dS_dT = similar(src)

	Threads.@threads for k in axes(src, 3)
		@inbounds for j in 1:nRho
			@inbounds for i in 1:nT
				dS, dT = if (i>1) &&(i<nT)
					log(src[i+1, j, k]) - log(src[i-1, j, k]),
					lnT[i+1] - lnT[i-1]
				elseif i==1
					log(src[i+1, j, k]) - log(src[i, j, k]),
					lnT[i+1] - lnT[i]
				elseif i==nT
					log(src[i, j, k]) - log(src[i-1, j, k]),
					lnT[i] - lnT[i-1]
				end
				dS_dT[i, j, k] = exp(log(src[i, j, k]) - lnT[i]) * dS / dT
			end
		end
	end

	tab.opacity[:dS_dT] = dS_dT

	dS_dT
end

function gradients!(eos::EquationOfState)
	@assert !is_internal_energy(eos)
	χₜ, χᵨ, cᵥ, dlnRoss_dlnT = _gradients(eos[:lnT], eos[:lnRho], eos[:lnEi], eos[:lnPg], eos[:lnRoss])

	eos.eos[:χₜ] = χₜ
	eos.eos[:χᵨ] = χᵨ
	eos.eos[:cᵥ] = cᵥ
	eos.eos[:dlnκ_dlnT] = dlnRoss_dlnT

	χₜ, χᵨ, cᵥ
end

function _gradients(lnT, lnRho, lnEi, lnPg, lnRoss)
	nT, nRho = length(lnT), length(lnRho)
	χₜ = similar(lnPg)
	χᵨ = similar(lnPg)
	cᵥ = similar(lnPg)
	dlnRoss_dlnT = similar(lnPg)
	dlnRoss_dlnRho = similar(lnPg)

	@inbounds for j in eachindex(lnRho)
		@inbounds for i in eachindex(lnT)
			dE, dT, dP, dRoss = if (i>1) &&(i<nT)
				lnEi[i+1, j] - lnEi[i-1, j],
				lnT[i+1    ] - lnT[i-1],
				lnPg[i+1, j] - lnPg[i-1, j],
				lnRoss[i+1, j] - lnRoss[i-1, j]
			elseif i==1
				lnEi[i+1, j] - lnEi[i, j],
				lnT[i+1]     - lnT[i],
				lnPg[i+1, j] - lnPg[i, j],
				lnRoss[i+1, j] - lnRoss[i, j]
			elseif i==nT
				lnEi[i, j] - lnEi[i-1, j],
				lnT[i]     - lnT[i-1],
				lnPg[i, j] - lnPg[i-1, j],
				lnRoss[i, j] - lnRoss[i-1, j]
			end
			cᵥ[i, j] = max(dE/dT * exp(lnEi[i, j] - lnT[i]), 1e-12)
			χₜ[i, j] = max(dP/dT, 1e-12)
			dlnRoss_dlnT[i, j] = dRoss / dT

			dR, dP, dRoss = if (j>1) &&(j<nRho)
				lnRho[j+1] - lnRho[j-1],
				lnPg[i, j+1] - lnPg[i, j-1],
				lnRoss[i, j+1] - lnRoss[i, j-1]
			elseif j==1
				lnRho[j+1] - lnRho[j],
				lnPg[i, j+1] - lnPg[i, j],
				lnRoss[i, j+1] - lnRoss[i, j]
			elseif j==nRho
				lnRho[j] - lnRho[j-1],
				lnPg[i, j] - lnPg[i, j-1],
				lnRoss[i, j] - lnRoss[i, j-1]
			end
			χᵨ[i, j] = max(dP/dR, 1e-12)
			dlnRoss_dlnRho[i, j] = dRoss / dR
			dlnRoss_dlnT[i, j] = dlnRoss_dlnT[i, j] - χₜ[i, j] / χᵨ[i, j]*dlnRoss_dlnRho[i, j]
		end
	end

	χₜ, χᵨ, cᵥ, dlnRoss_dlnT
end

# ============================================================================
# Thermodynamic quantities for pre-computation
# ============================================================================

function add_thermodynamics!(eos::EquationOfState)
	χₜ, χᵨ, cᵥ = gradients!(eos)
	Q, cₚ, ∇ₐ, μ, S = _thermodynamics(χₜ, χᵨ, cᵥ, eos[:lnT], eos[:lnRho], eos[:lnEi], eos[:lnPg])
	eos.eos[:Q] = Q
	eos.eos[:cₚ] = cₚ
	eos.eos[:∇ₐ] = ∇ₐ
	eos.eos[:μ] = μ
	eos.eos[:S] = S

	nothing
end

function _thermodynamics(χₜ, χᵨ, cᵥ, lnT, lnRho, lnEi, lnPg)
	Q = χₜ ./ χᵨ

	cₚ = similar(cᵥ)
	∇ₐ = similar(cᵥ)
	μ = similar(cᵥ)
	S = zeros(size(cᵥ))
	for j in eachindex(lnRho)
		for i in eachindex(lnT)
			P = exp(lnPg[i, j])
			ρ = exp(lnRho[j])
			T = exp(lnT[i])
			cₚ[i, j] = min(1e12, max(1e-12, cᵥ[i, j] + P/(ρ*T) * χₜ[i, j]^2 / χᵨ[i, j]))
			∇ₐ[i, j] = min(1, max(1e-12, P/(ρ*T) .* Q[i, j] ./ cₚ[i, j]))
			μ[i, j] = min(1e12, max(1e-12, ρ * KBoltzmann * T / (P * m_u)))
		end
	end

	S[1, 1] = ideal_entropy(exp(lnT[1]), exp(lnRho[1]), μ[1, 1])
	for i in 2:size(S, 1)
		dE = lnEi[i, 1] - lnEi[i-1, 1]
		T = (exp(lnT[i]) +  exp(lnT[i-1])) / 2
		E = (exp(lnEi[i, 1]) +  exp(lnEi[i-1, 1])) / 2

		dS = E / T * dE
		S[i, 1] = S[i-1, 1] + dS
	end

	T = exp.(lnT)
	for j in 2:size(S, 2)
		dE = lnEi[:, j] - lnEi[:, j-1]
		dRho = lnRho[j] - lnRho[j-1]

		rho = (exp(lnRho[j]) + exp(lnRho[j-1])) / 2
		P = (exp.(lnPg[:, j]) + exp.(lnPg[:, j-1])) / 2
		E = (exp.(lnEi[:, j]) + exp.(lnEi[:, j-1])) / 2

		t1 = E ./ T .* dE
		t2 = P ./ (rho .* T) .* dRho

		dS = t1 - t2
		S[:, j] = S[:, j-1] .+ dS
	end

	Q, cₚ, ∇ₐ, μ, S
end
