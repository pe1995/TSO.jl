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

opacity(opa::ExtendedOpacity) = opa.opa
wavelength(opa::ExtendedOpacity) = wavelength(opa.opa)

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

lookup(eos::ExtendedEoS, opa::OpacityTable, args...; kwargs...) = lookup(eos.eos, opa, args...; kwargs...)
lookup(eos::ExtendedEoS, args...; kwargs...) = extended_lookup(eos, args...; kwargs...)


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
lookup(eos, opa::ExtendedOpacity, args...; kwargs...) = extended_lookup(eos, opa, args...; kwargs...)








function ideal_entropy(T, rho, mu)    
    m_particle = mu * m_u
    R_spec = KBoltzmann / m_particle
    
    factor = (2.0 * pi * m_particle * KBoltzmann * T) / (HPlanck^2)
    n_quantum = factor^1.5
    n = rho / m_particle
    
    # Sackur-Tetrode Equation
    R_spec * (log(n_quantum / n) + 2.5)
end







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