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
    @assert is_internal_energy(eos)
    eemin,ee_max,_,_ = exp.(limits(eos))

    ee_min = if isnothing(ee_min)
        eemin
    else
        ee_min
    end

    # initial point
    t0  = exp(first(star_point.lnT))
    ρ0  = exp(first(star_point.lnρ))
    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
    g0  = exp10(star_point.logg)
    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))
    z0  = first(star_point.z)

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
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
        
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
            p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
        
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
    model = deepcopy(model_in)

    # First make sure the z scale is oriented in the correct direction
    if first(model.lnT) < last(model.lnT)
        @warn "Flipping z axis..."
        flip!(model)
    end

    if first(model.z) > last(model.z)
        model.z .= - model.z
    end

    start_point = pick_point(model, 1)
    end_point = pick_point(model, length(model.z))
    
    adiabat(start_point, end_point, eos; kwargs...)
end