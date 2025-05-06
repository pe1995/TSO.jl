"""
    adiabat(model; kwargs...)

Construct a 1D adiabat just from start and end points.
"""
function adiabat(star_point, end_point, eos_in::SqEoS; 
                    nz=440, n_iter=30, dlnd=0.05, padding=0.05, ee_min=nothing)
    eos = @axed eos_in
    
    zbottom = first(star_point.z)
    ztop = first(end_point.z)

    # Find the initial energy from the EoS. Assert that the EoS is in the 
    # correct units for this
    ee_min, ee_max = if is_internal_energy(eos)
        eemin,ee_max,_,_ = exp.(limits(eos))

        ee_min = if isnothing(ee_min)
            eemin
        else
            ee_min
        end

        ee_min, ee_max
    else
        eemin,ee_max = exp(minimum(eos.eos.lnEi)), exp(maximum(eos.eos.lnEi))

        ee_min = if isnothing(ee_min)
            eemin
        else
            ee_min
        end

        ee_min, ee_max
    end

    # initial point
    t0 = exp(first(star_point.lnT))
    ρ0 = exp(first(star_point.lnρ))
    g0 = exp10(star_point.logg)
    z0  = first(star_point.z)

    ee0, p0 = if is_internal_energy(eos)
        exp(lookup(eos, :lnEi, log(ρ0), log(t0))), exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))
    else
        exp(lookup(eos, :lnEi, log(ρ0), log(t0))), exp(lookup(eos, :lnPg, log(ρ0), log(t0)))
    end

    # storage arrays
    z  = zeros(2*nz)
    d  = zeros(2*nz)
    ee = zeros(2*nz)
    t  = zeros(2*nz)
    p  = zeros(2*nz)

    # integrate from the starting point upwards, until z is larger than ztop
    # or the maximum number of points is reached
    dee    = 0
    z[nz]  = z0
    d[nz]  = ρ0  
    ee[nz] = ee0 
    t[nz]  = t0  
    p[nz]  = p0  
    i_end = 2*nz

    δz_sponge = abs(first(end_point.z) - first(star_point.z)) * padding
    
    # upwards
    for i in nz+1:2*nz
        ee[i] = max(ee[i-1]-dee,ee_min)
        
        for iter in 1:n_iter
            d[i] = d[i-1]*exp(-dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = if is_internal_energy(eos)
                exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
            else
                exp(lookup(eos, :lnPg, log(d[i]), log(t[i])))
            end
        
            dp = p[i-1] - p[i]
            da = (d[i-1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i-1] + dz
            pd    = (p[i-1] + p[i]) / (d[i-1] + d[i])
            dee   = pd*dlnd
            ee[i] = max(ee[i-1]-dee, ee_min)
        end

        if z[i] > ztop + δz_sponge
            i_end = i
            break
        end
    end

    # downwards
    i_start=1
    for i in nz:-1:1
        ee[i] = min(ee[i+1]+dee, ee_max)
        
        for iter in 1:n_iter
            d[i] = d[i+1]*exp(dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = if is_internal_energy(eos)
                exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
            else
                exp(lookup(eos, :lnPg, log(d[i]), log(t[i])))
            end        
            dp = p[i] - p[i+1]
            da = (d[i+1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i+1] - dz
            pd    = (p[i+1] + p[i]) / (d[i+1] + d[i])
            dee   = pd*dlnd
            ee[i] = min(ee[i+1]+dee, ee_max)
        end

        if z[i] < zbottom - δz_sponge
            i_start = i
            break
        end
    end
    
    z = z[i_start:i_end]
    d = d[i_start:i_end]
    ee = ee[i_start:i_end]
    t = t[i_start:i_end] 
    p = p[i_start:i_end]  

    Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=star_point.logg)
end

"""
    adiabat(model; kwargs...)

Construct a 1D adiabatic model by integrating from the bottom to the top. 
A new height scale is constructed automatically. The integration is 
executed for max. nz points in both directions. When the bottom
of the average3D model is reached, it stops automatically.
"""
function adiabat(model_in, eos::SqEoS; kwargs...)
    # First make sure the z scale is oriented in the correct direction
    model = flip(model_in)

    start_point = pick_point(model, 1)
    end_point = pick_point(model, length(model))
    
    adiabat(start_point, end_point, eos; kwargs...)
end







#= Adiabatic extrapolation =#

function adiabatDown(start_point, eos_in::SqEoS; kwargs...)
    # initial point
    t0  = exp(first(start_point.lnT))
    ρ0  = exp(first(start_point.lnρ))
    g0  = exp10(start_point.logg)
    z0  = first(start_point.z)

    adiabatDown(t0, ρ0, z0, start_point.logg, eos_in; kwargs...)
end

adiabatDown(t0, ρ0, z0, logg, eos_in::SqEoS; kwargs...) = is_internal_energy(@axed(eos_in)) ? 
    adiabatDown_ee(t0, ρ0, z0, logg, eos_in; kwargs...) : 
    adiabatDown_t(t0, ρ0, z0, logg, eos_in; kwargs...)

function adiabatDown_ee(t0, ρ0, z0, logg, eos_in::SqEoS; 
                    nz=10000, n_iter=30, dlnd=0.075, ee_min=nothing, z_end=-Inf)
    eos = @axed eos_in

    # Find the initial energy from the EoS. Assert that the EoS is in the 
    # correct units for this
    @assert is_internal_energy(eos)
    eemin,ee_max,_,_ = exp.(limits(eos))

    ee_min = if isnothing(ee_min)
        eemin
    else
        ee_min
    end

    # initial point
    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
    g0  = exp10(logg)
    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))

    # storage arrays
    z  = zeros(nz)
    d  = zeros(nz)
    ee = zeros(nz)
    t  = zeros(nz)
    p  = zeros(nz)

    # integrate from the starting point upwards, until z is larger than ztop
    # or the maximum number of points is reached
    dee    = 0
    z[nz]  = z0
    d[nz]  = ρ0  
    ee[nz] = ee0 
    t[nz]  = t0  
    p[nz]  = p0  

    # downwards
    i_start=1
    for i in nz-1:-1:1
        ee[i] = min(ee[i+1]+dee, ee_max)
        
        for iter in 1:n_iter
            d[i] = d[i+1]*exp(dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
        
            dp = p[i] - p[i+1]
            da = (d[i+1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i+1] - dz
            pd    = (p[i+1] + p[i]) / (d[i+1] + d[i])
            dee   = pd*dlnd
            ee[i] = min(ee[i+1]+dee, ee_max)
        end

        # check if there is an end 
        if z[i] <= z_end
            i_start = i
            break
        end
    end

    z = z[i_start:end]
    d = d[i_start:end]
    t = t[i_start:end]
    ee = ee[i_start:end]
    
    Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=logg)
end

function adiabatDown_t(t0, ρ0, z0, logg, eos_in::SqEoS; 
                    nz=10000, n_iter=30, dlnd=0.075, ee_min=nothing, z_end=-Inf)
    eos = @axed eos_in

    # Find the initial energy from the EoS. Assert that the EoS is in the 
    # correct units for this
    @assert !is_internal_energy(eos)
    eemin,ee_max = exp(minimum(eos.eos.lnEi)), exp(maximum(eos.eos.lnEi))

    ee_min = if isnothing(ee_min)
        eemin
    else
        ee_min
    end

    # initial point
    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
    g0  = exp10(logg)
    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(t0)))

    # storage arrays
    z  = zeros(nz)
    d  = zeros(nz)
    ee = zeros(nz)
    t  = zeros(nz)
    p  = zeros(nz)

    # integrate from the starting point upwards, until z is larger than ztop
    # or the maximum number of points is reached
    dee    = 0
    z[nz]  = z0
    d[nz]  = ρ0  
    ee[nz] = ee0 
    t[nz]  = t0  
    p[nz]  = p0  

    # downwards
    i_start=1
    for i in nz-1:-1:1
        ee[i] = min(ee[i+1]+dee, ee_max)
        
        for iter in 1:n_iter
            d[i] = d[i+1]*exp(dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(t[i])))
        
            dp = p[i] - p[i+1]
            da = (d[i+1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i+1] - dz
            pd    = (p[i+1] + p[i]) / (d[i+1] + d[i])
            dee   = pd*dlnd
            ee[i] = min(ee[i+1]+dee, ee_max)
        end

        # check if there is an end 
        if z[i] <= z_end
            i_start = i
            break
        end
    end

    z = z[i_start:end]
    d = d[i_start:end]
    t = t[i_start:end]
    ee = ee[i_start:end]
    
    Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=logg)
end




function adiabatUp(start_point, eos_in::SqEoS; kwargs...)
    # initial point
    t0  = exp(first(start_point.lnT))
    ρ0  = exp(first(start_point.lnρ))
    g0  = exp10(start_point.logg)
    z0  = first(start_point.z)

    adiabatUp(t0, ρ0, z0, start_point.logg, eos_in; kwargs...)
end

adiabatUp(t0, ρ0, z0, logg, eos_in::SqEoS; kwargs...) = is_internal_energy(@axed(eos_in)) ? 
    adiabatUp_ee(t0, ρ0, z0, logg, eos_in; kwargs...) : 
    adiabatUp_t(t0, ρ0, z0, logg, eos_in; kwargs...)

function adiabatUp_ee(t0, ρ0, z0, logg, eos_in::SqEoS; 
                    nz=10000, n_iter=30, dlnd=0.075, ee_min=nothing, z_end=Inf)
    eos = @axed eos_in

    # Find the initial energy from the EoS. Assert that the EoS is in the 
    # correct units for this
    @assert is_internal_energy(eos)
    eemin,ee_max,_,_ = exp.(limits(eos))

    ee_min = if isnothing(ee_min)
        eemin
    else
        ee_min
    end

    # initial point
    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
    g0  = exp10(logg)
    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))

    # storage arrays
    z  = zeros(nz)
    d  = zeros(nz)
    ee = zeros(nz)
    t  = zeros(nz)
    p  = zeros(nz)

    # integrate from the starting point upwards, until z is larger than ztop
    # or the maximum number of points is reached
    dee    = 0
    z[1]  = z0
    d[1]  = ρ0  
    ee[1] = ee0 
    t[1]  = t0  
    p[1]  = p0
    i_end = nz


    # upwards
    for i in 2:nz
        ee[i] = max(ee[i-1]-dee,ee_min)
        
        for iter in 1:n_iter
            d[i] = d[i-1]*exp(-dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
        
            dp = p[i-1] - p[i]
            da = (d[i-1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i-1] + dz
            pd    = (p[i-1] + p[i]) / (d[i-1] + d[i])
            dee   = pd*dlnd
            ee[i] = max(ee[i-1]-dee, ee_min)
        end

        if z[i] > z_end
            i_end = i
            break
        end
    end

    z = z[1:i_end]
    d = d[1:i_end]
    t = t[1:i_end]
    ee = ee[1:i_end]
    
    Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=logg)
end

function adiabatUp_t(t0, ρ0, z0, logg, eos_in::SqEoS; 
                    nz=10000, n_iter=30, dlnd=0.075, ee_min=nothing, z_end=Inf)
    eos = @axed eos_in

    # Find the initial energy from the EoS. Assert that the EoS is in the 
    # correct units for this
    @assert !is_internal_energy(eos)
    eemin,ee_max = exp(minimum(eos.eos.lnEi)), exp(maximum(eos.eos.lnEi))

    ee_min = if isnothing(ee_min)
        eemin
    else
        ee_min
    end

    # initial point
    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
    g0  = exp10(logg)
    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(t0)))

    # storage arrays
    z  = zeros(nz)
    d  = zeros(nz)
    ee = zeros(nz)
    t  = zeros(nz)
    p  = zeros(nz)

    # integrate from the starting point upwards, until z is larger than ztop
    # or the maximum number of points is reached
    dee    = 0
    z[1]  = z0
    d[1]  = ρ0  
    ee[1] = ee0 
    t[1]  = t0  
    p[1]  = p0  
    i_end = nz

    # upwards
    for i in 2:nz
        ee[i] = max(ee[i-1]-dee,ee_min)
        
        for iter in 1:n_iter
            d[i] = d[i-1]*exp(-dlnd)
            t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(t[i])))
        
            dp = p[i-1] - p[i]
            da = (d[i-1] + d[i]) /2.
            dz = dp / (g0*da)
            
            z[i]  = z[i-1] + dz
            pd    = (p[i-1] + p[i]) / (d[i-1] + d[i])
            dee   = pd*dlnd
            ee[i] = max(ee[i-1]-dee, ee_min)
        end

        if z[i] > z_end
            i_end = i
            break
        end
    end

    z = z[1:i_end]
    d = d[1:i_end]
    t = t[1:i_end]
    ee = ee[1:i_end]
    
    Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=logg)
end





stitch(model, adiabat) = begin
	z = [adiabat.z..., model.z...]
	mask = unique(i -> z[i], eachindex(z))
	Model1D(
		z=[adiabat.z..., model.z...][mask], 
		lnρ=[adiabat.lnρ..., model.lnρ...][mask], 
		lnT=[adiabat.lnT..., model.lnT...][mask], 
		lnEi=[adiabat.lnEi..., model.lnEi...][mask], 
		logg=model.logg
    )
end



function adiabatic_extrapolation(model, eos, dz; kwargs...)
    model_flip = flip(model)
    i_top = pick_point(model_flip, 1)

    z_end = first(model_flip.z) - dz
    a = adiabatDown(i_top, eos, z_end=z_end; kwargs...)
    ad_opt = @optical a eos
    
	flip!(ad_opt, depth=false)
	model_flip.lnEi .= lookup(eos, :lnEi, model_flip.lnρ, model_flip.lnT)
	
    model = flip(
        Model1D(
            z=[ad_opt.z..., model_flip.z...], 
            lnρ=[ad_opt.lnρ..., model_flip.lnρ...], 
            lnT=[ad_opt.lnT..., model_flip.lnT...], 
            lnEi=[ad_opt.lnEi..., model_flip.lnEi...], 
            logg=model_flip.logg
        ), 
        depth=true
    )
	uniform_z = range(minimum(model.z), maximum(model.z), length=length(model.z)) 
	flip(
		interpolate_to(model, in_log=false, z=uniform_z |> collect),
		depth=true
	)
end

function adiabatic_extrapolation(model, eos; τ_extrapolate=1.0, τ_target=8.0, iter_max=500, kwargs...)
    model_flip = @optical flip(model) eos

	# crop the model until τ=τ_extrapolate
	model_opt_cut = flip(
        cut(
            flip(model_flip, depth=true), 
            τ=[minimum(model_flip.τ), exp10(τ_extrapolate)]
        ), 
        depth=false
    )
	
	tau_bottom = Ref(-Inf)
	mloc = deepcopy(model_opt_cut)
	ic = 0
	while (tau_bottom[] < τ_target) && (ic < iter_max)
		ic += 1

		i_top = pick_point(mloc, 1)
		a = adiabatDown(i_top, eos; nz=10, dlnd=0.01, kwargs...)
		mloc = stitch(mloc, a)

		τ = rosseland_optical_depth(eos, flip(mloc, depth=true))
		tau_bottom[] = maximum(log10.(τ)) 
	end
	
	mloc.lnEi .= lookup(eos, :lnEi, mloc.lnρ, mloc.lnT)
	uniform_z = range(minimum(mloc.z), maximum(mloc.z), length=length(mloc.z)) 
	@optical flip(
		interpolate_to(mloc, in_log=false, z=uniform_z |> collect),
		depth=true
	) eos
end

function reverse_adiabatic_extrapolation(model, rho_bottom, T_bottom, eos; τbottom=1.0, iter_max=500, kwargs...)
    model_flip = flip(model)

	# crop the model until τ=τbottom
	model_opt_cut = flip(cut(flip(model_flip, depth=true), τ=[minimum(model_flip.τ), exp10(τbottom)]), depth=false)
	
	mloc = deepcopy(model_opt_cut)

	tau_top = Ref(Inf)
	ic = 0
	a = Ref{Any}(0)
	while (tau_top[] > τbottom) && (ic < iter_max)
		ic += 1
		a[] = adiabatUp(
			T_bottom, rho_bottom, 0.0, model.logg, eos; z_end=Inf, nz=100
		)
		τ = rosseland_optical_depth(eos, flip(a[], depth=true))
		tau_top[] = minimum(log10.(τ)) 
	end
	ao = @optical a[] eos
	a = flip(
		cut(
			flip(ao, depth=true), 
			lnρ=[maximum(model_opt_cut.lnρ), maximum(ao.lnρ)]
		), 
		depth=false
	)

	ref_z = mloc.z[1]
	a.z .-= a.z[end] 
	a.z .+= ref_z
	mloc = stitch(mloc, a)
	
	mloc.lnEi .= lookup(eos, :lnEi, mloc.lnρ, mloc.lnT)
	uniform_z = range(minimum(mloc.z), maximum(mloc.z), length=length(mloc.z)) 
	o = @optical flip(
		interpolate_to(mloc, in_log=false, z=uniform_z |> collect),
		depth=true
	) eos

	optical_surface!(o)

	o
end

