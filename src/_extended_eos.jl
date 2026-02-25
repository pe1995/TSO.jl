# ============================================================================
# Extended EoS -- Flexible extension of the EoS for arbitrary quantities
# ============================================================================
@kwdef struct ExtendedEoS
	eos 
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
extended(eos::RegularEoSTable) = ExtendedEoS(eos=eos)

get_data(eos::ExtendedEoS, opa::ExtendedOpacity, var::Symbol) = begin
    Typ = eltype(eos.eos.lnPg)
    if var in fieldnames(typeof(eos.eos))
        return getfield(eos.eos, var)::AbstractArray{Typ}
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

function extended_lookup(eos::ExtendedEoS, what, pressure_var, energy_var; maxiter=200, tol=1e-8)
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
	elseif what in keys(eos.extensions)
		x = eos.eos.lnT
		y = eos.eos.lnRho
		z = eos.extensions[what]
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
add_log_interpolated(var::Symbol) = push!(LOG_INTERPOLATED_VARS, var)

weights(eos::ExtendedEoS, lnrho::AbstractArray, lnT::AbstractArray) = begin
    grid_T = eos.eos.lnT
    grid_Rho = eos.eos.lnRho
    
    #T = eltype(grid_T)
    T = eltype(lnT)
    
    coefs_T = Array{InterpCoefs{T}}(undef, size(lnT))
    coefs_Rho = Array{InterpCoefs{T}}(undef, size(lnrho))
    
    @inbounds for j in eachindex(lnT)
        coefs_T[j] = linear_interpolation_weights(grid_T, lnT[j])
        coefs_Rho[j] = linear_interpolation_weights(grid_Rho, lnrho[j])
    end
    
    coefs_Rho, coefs_T
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
sample(eos::ExtendedEoS, opa::RegularOpacityTable, args...; kwargs...) =  sample(eos, extended(opa), args...; kwargs...)
sample(eos::ExtendedEoS, opa::ExtendedOpacity, variables, lnrho::Number, lnT::Number, args...) = begin
    res = sample(eos, opa, variables, [lnrho], [lnT], args...)
    return map(r -> ndims(r) == 1 ? r[1] : vec(r[:, 1]), res)
end

function sample(eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, lnrho::AbstractArray{T}, lnT::AbstractArray{T}, args...) where {T<:AbstractFloat}
    results = map(variables) do var
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
    sample!(results, eos, opa, variables, lnrho, lnT, args...)
    results
end

sample!(results::Union{Tuple, AbstractArray}, eos::RegularEoSTable, args...; kwargs...) = sample!(results, extended(eos), args...; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::RegularOpacityTable, args...; kwargs...) = sample!(results, eos, extended(opa), args...; kwargs...)
sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, lnrho::AbstractArray{T}, lnT::AbstractArray{T}, args...) where {T<:AbstractFloat} = begin
    coefs_Rho, coefs_T = weights(eos, lnrho, lnT)
    sample!(results, eos, opa, variables, coefs_Rho, coefs_T, args...)
    results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        data = get_data(eos, opa, var)
        is_log = _is_log_interpolated(var)
        _interpolate_var!(results[i], data, coefs_Rho, coefs_T, is_log)
    end
    
    return results
end

function sample!(results::Union{Tuple, AbstractArray}, eos::ExtendedEoS, opa::ExtendedOpacity, variables::Union{Tuple, AbstractArray}, coefs_Rho::AbstractArray{InterpCoefs{T2}}, coefs_T::AbstractArray{InterpCoefs{T2}}, k::Union{Integer, AbstractUnitRange}) where {T2<:AbstractFloat}
    for (i, var) in enumerate(variables)
        data = get_data(eos, opa, var)
        is_log = _is_log_interpolated(var)
        _interpolate_var!(results[i], data, coefs_Rho, coefs_T, k, is_log)
    end
    
    return results
end

# ============================================================================
# Core interpolation
# ============================================================================

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
            for k in 1:n3
                result[k + (j-1)*n3] = _interpolate_point_3d(data, cr, ct, k, is_log)
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
                result[j] = _interpolate_point_3d(data, coefs_Rho[j], coefs_T[j], k, is_log)
            end
        else
            nk = length(k)
            @inbounds for j in 1:length(coefs_T)
                cr = coefs_Rho[j]
                ct = coefs_T[j]
                for (ik, kidx) in enumerate(k)
                    result[ik + (j-1)*nk] = _interpolate_point_3d(data, cr, ct, kidx, is_log)
                end
            end
        end
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

@inline function _interpolate_point_3d(data::AbstractArray{T}, cr::InterpCoefs{T2}, ct::InterpCoefs{T2}, k::Int, is_log::Bool) where {T<:AbstractFloat, T2<:AbstractFloat}
    it, ir = ct.idx, cr.idx
    w00 = ct.w_low  * cr.w_low
    w10 = ct.w_high * cr.w_low
    w01 = ct.w_low  * cr.w_high
    w11 = ct.w_high * cr.w_high
    
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
	eos = eos_extended.eos
	aos = @axed eos
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
	eos = eos_extended.eos
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