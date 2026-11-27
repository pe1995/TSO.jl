# ============================================================================
# Interpolate tables
# ============================================================================

_float_fields(d) = isnothing(d) ? Symbol[] : sort([k for k in keys(d) if eltype(d[k]) <: AbstractFloat])
_eos_fields(tab::EquationOfState) = _float_fields(tab.eos)

"""
    switch_energy(tab; upsample=-1, newE=nothing)

Switch the energy axis of `tab` (temperature <-> internal energy). All quantities
are interpolated column by column, opacities in log.
"""
function switch_energy(tab::EquationOfState; upsample=-1, newE=nothing)
    newname = dependent_energy_variable(tab)
    oldname = energy_variable(tab)
    newE = isnothing(newE) ? pick_axis(tab; upsample=upsample) : newE

    nE, nρ = length(newE), tab.density_axis.length
    nOld = tab.energy_axis.length
    T = eltype(tab.density_axis.values)
    newAxisOldE = tab[newname]

    eos_in = Dict{Symbol, AbstractArray}(f => tab.eos[f] for f in _eos_fields(tab) if f != newname)
    oldaxis = similar(newAxisOldE)
    puff_up!(oldaxis, tab.energy_axis.values, 1)
    eos_in[oldname] = oldaxis
    eos_out = Dict{Symbol, AbstractArray}(f => zeros(T, nE, nρ) for f in keys(eos_in))

    parts = [p for p in (:opacity, :scattering) if !isnothing(getfield(tab, p))]
    opa_out = Dict(p => Dict{Symbol, AbstractArray}(f => zeros(T, nE, nρ, size(getfield(tab, p)[f])[3:end]...) for f in _float_fields(getfield(tab, p))) for p in parts)

    x = zeros(nOld)
    y = zeros(nOld)
    mask = zeros(Int, nOld)

    @inbounds for i in 1:nρ
        mask .= sortperm(view(newAxisOldE, :, i))
        x .= Interpolations.deduplicate_knots!(newAxisOldE[mask, i], move_knots=true)

        for (f, v) in eos_in
            y .= @view(v[mask, i])
            ip = linear_interpolation(x, y, extrapolation_bc=Line())
            eos_out[f][:, i] .= ip.(newE)
        end

        for p in parts, (f, out) in opa_out[p]
            v = getfield(tab, p)[f]
            for j in axes(out, 3)
                y .= log.(@view(v[mask, i, j]))
                ip = linear_interpolation(x, y, extrapolation_bc=Line())
                out[:, i, j] .= ip.(newE)
            end
        end
    end

    _new_opacity(p) = begin
        o = getfield(tab, p)
        OpacityData(Dict{Symbol, AbstractArray}(f => exp.(v) for (f, v) in opa_out[p]), deepcopy(o.λ), deepcopy(o.weights), o.binned, o.optical_depth)
    end

    EquationOfState(deepcopy(tab.density_axis.values), newname => newE;
        eos=eos_out,
        opacity=has_opacity(tab) ? _new_opacity(:opacity) : nothing,
        scattering=has_scattering(tab) ? _new_opacity(:scattering) : nothing,
        meta=tab.meta
    )
end

"""
    interpolate_2D(tab; lnRho, lnT/lnEi)

Interpolate (and extrapolate) all quantities of `tab` to new values of its axes.
"""
function interpolate_2D(tab::EquationOfState; lnRho, newE...)
    name = only(keys(newE))
    name == energy_variable(tab) || error("interpolate_2D only changes the values of the axes, $(name) is not the energy axis of the table.")
    newAxis = newE[name]
    ee, rr = meshgrid(newAxis, lnRho)

    resample(f) = sample(tab, (f,), rr, ee; extrapolate=true)[1]
    eos = Dict{Symbol, AbstractArray}(f => resample(f) for f in _eos_fields(tab))

    _new_opacity(p) = begin
        o = getfield(tab, p)
        t = EquationOfState(tab; opacity=o)
        OpacityData(Dict{Symbol, AbstractArray}(f => permutedims(sample(t, (f,), rr, ee; extrapolate=true)[1], ndims(o[f]) == 3 ? (2, 3, 1) : (1, 2)) for f in _float_fields(o)),
                    deepcopy(o.λ), deepcopy(o.weights), o.binned, o.optical_depth)
    end

    new = EquationOfState(collect(lnRho), name => collect(newAxis);
        eos=eos,
        opacity=has_opacity(tab) ? _new_opacity(:opacity) : nothing,
        scattering=has_scattering(tab) ? _new_opacity(:scattering) : nothing,
        meta=tab.meta
    )
    fill_nan!(new)
    new
end

pick_axis(tab::EquationOfState; upsample=-1) = begin
    values = tab[dependent_energy_variable(tab)]
    sort(collect(range(minimum(values), maximum(values), length=upsample == -1 ? tab.energy_axis.length : upsample)))
end

function puff_up!(Tnew::AbstractArray{T,2}, Told::Vector, dim=1) where {T}
    @inbounds for j in axes(Tnew, 2)
        @inbounds for i in axes(Tnew, 1)
            Tnew[i, j] = dim==1 ? Told[i] : Told[j]
        end
    end
    Tnew
end

# ============================================================================
# NaN handling
# ============================================================================

@inline nan_or_inf(a) = isnan(a) | !isfinite(a)
@inline interpolate_at(x, y, x0; bc=Flat()) = linear_interpolation(Interpolations.deduplicate_knots!(deepcopy(x), move_knots=true), y, extrapolation_bc=bc).(x0)
@inline function _lerp_flat(x_knots, y_knots, n, x_query)
    @inbounds x_query <= x_knots[1] && return y_knots[1]
    @inbounds x_query >= x_knots[n] && return y_knots[n]

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

function _interpolate_at!(y_slice, x_coords, row_mask, buf_x, buf_y;
                           is_log=false, fallback=1.0f-30)
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
        @inbounds for i in eachindex(y_slice)
            if row_mask[i]
                y_slice[i] = fallback
            end
        end
    end
    return true
end

"""
    nan_mask(tab; eos_fields=all, opacities=true)

Cells where any of the given quantities is NaN, infinite (or not positive for opacities).
"""
function nan_mask(tab::EquationOfState; eos_fields=_eos_fields(tab), opacities=true)
    mask = falses(size(tab)...)
    for f in eos_fields
        mask .|= nan_or_inf.(tab.eos[f])
    end
    if opacities
        for p in (:opacity, :scattering), o in (getfield(tab, p),)
            isnothing(o) && continue
            for f in _float_fields(o)
                bad = x -> nan_or_inf(x) || x <= 0
                v = o[f]
                mask .|= ndims(v) == 2 ? bad.(v) : dropdims(any(bad, v, dims=3), dims=3)
            end
        end
    end
    Base.convert(Array{Bool, 2}, mask)
end

"""
    fill_nan!(tab; eos_fields=all, opacities=true, fallback=1f-30)

Replace NaN and Inf in the table by interpolating along the energy axis, then
along the density axis. Returns the mask of the cells that were replaced.
"""
function fill_nan!(tab::EquationOfState; eos_fields=_eos_fields(tab), opacities=true, fallback=1.0f-30)
    x_rho = tab.density_axis.values
    x_e   = tab.energy_axis.values

    nan_mask = falses(length(x_e), length(x_rho))
    lfallback = log(fallback)

    eos_arrays = [tab.eos[f] for f in eos_fields]
    opa_arrays = AbstractArray[]
    if opacities
        for p in (:opacity, :scattering)
            o = getfield(tab, p)
            isnothing(o) && continue
            append!(opa_arrays, [o[f] for f in _float_fields(o)])
        end
    end

    for v in eos_arrays
        _fill_nan!(nan_mask, v, x_e, x_rho, false, lfallback)
    end
    for v in opa_arrays
        _fill_nan!(nan_mask, v, x_e, x_rho, true, fallback)
    end

    return Base.convert(Array{Bool, 2}, nan_mask)
end

function _fill_nan!(nan_mask, v, x_e, x_rho, is_log, fallback)
    col_mask = Vector{Bool}(undef, length(x_e))
    buf_xe   = similar(x_e)
    buf_ye   = similar(x_e)
    row_mask = Vector{Bool}(undef, length(x_rho))
    buf_x    = similar(x_rho)
    buf_y    = similar(x_rho)

    @inbounds for k in CartesianIndices(size(v)[3:end])
        for j in axes(v, 2)
            had = _interpolate_at!(view(v, :, j, k), x_e, col_mask, buf_xe, buf_ye;
                                   is_log=is_log, fallback=fallback)
            had && (view(nan_mask, :, j) .|= col_mask)
        end
        for i in axes(v, 1)
            _interpolate_at!(view(v, i, :, k), x_rho, row_mask, buf_x, buf_y;
                             is_log=is_log, fallback=fallback)
        end
    end
    v
end

"""
    set_limits!(tab; small=1e-30, large=1e30)

Clamp all quantities of the table to [small, large] (the EoS in log), NaN is set to `small`.
"""
function set_limits!(tab::EquationOfState; small=1e-30, large=1e30)
    lsmall, llarge = log(small), log(large)
    for f in _eos_fields(tab)
        _set_limits!(tab.eos[f], lsmall, llarge)
    end
    for p in (:opacity, :scattering)
        o = getfield(tab, p)
        isnothing(o) && continue
        for f in _float_fields(o)
            _set_limits!(o[f], small, large)
        end
    end
    tab
end

function _set_limits!(v, small, large)
    @inbounds for i in eachindex(v)
        v[i] = check(v[i], small) ? v[i] : small
        v[i] = check_large(v[i], large) ? v[i] : large
    end
    v
end

"""
    smoothAccumulate!(tab; spline=true)

Make the internal energy monotonic in temperature (on the temperature grid).
"""
smoothAccumulate!(tab::EquationOfState; spline=true) = begin
    @assert !is_internal_energy(tab)
    _smooth_accumulate!(tab[:lnEi], tab[:lnT], spline)
    tab
end

function _smooth_accumulate!(lnEi, lnT, spline)
    x = copy(lnT)
	y = similar(lnT)

    for i in axes(lnEi, 2)
        y .= reverse(lnEi[:, i])
        y .= reverse(accumulate(min, y))

        lnEi[:, i] .= if spline
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
end

check(v, small) = !isnan(v) & (v >= small)
check_small = check
check_large(v, large) = !isnan(v) & (v < large)

