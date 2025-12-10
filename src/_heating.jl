#================================ Solvers ====================================#

struct AngleHeatingSolver{T<:AbstractFloat} <:RadiativeTransferSolver
    μ  ::T
    ω  ::T
    z  ::Vector{T}
    χ  ::Vector{T}
    S  ::Vector{T}
    Q  ::Array{T}
    
    # Process quantities
    dtau1 ::Array{T}
    dtau2 ::Array{T}
    ex0 ::Array{T}
    ex1 ::Array{T}
    ex2 ::Array{T}
    dsdtau1 ::Array{T}
    dsdtau2 ::Array{T}
    d1sdtau1 ::Array{T}
    d2sdtau2 ::Array{T}
end

struct HeatingSolver{T<:AbstractFloat, A<:AxedEoS, O<:OpacityTable} <:RadiativeTransferSolver
    eos          ::A
    opacities    ::O
    angles       ::Vector{T}
    weights      ::Vector{T}
    bins         ::Int
    model        ::AbstractModel
    χ            ::Vector{T}
    S            ::Vector{T}
    Q            ::Array{T}
    F            ::Array{T}
    μ_solver     ::Vector{AngleHeatingSolver{T}}
end

struct HeatingField{A<:AbstractModel, T<:AbstractFloat, N} <:AbstractRadiationField
    model::A
    values::Array{T, N}
    flux::Array{T, N}
end

struct InterpolatedHeating{R<:RadiationField, I} <:AbstractRadiationField
    radiation::R
    interpolation_function::Vector{I}
end






Heating(
    model::AbstractModel; 
    eos::AxedEoS, 
    opacities::AbstractBinnedOpacities, 
    angles=sort(labatto_4angles), 
    weights=J_weights(labatto_4angles)
) = begin
    opacities = opacities.opacities
    bins = size(opacities.κ, 3)

    T2 = Float64 

    # Make sure the model has the correct orientation
    m = convert_model(model, T2)
    flip!(m)

    mask = sortperm(angles, rev=true)
    angles = Base.convert.(T2, angles[mask])
    weights = Base.convert.(T2, weights[mask])

    # Create an angular solver for every angle in angles, with the correct weight
    angle_solver = [
        AngleHeatingSolver(
            angles[i], weights[i], 
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z)),
            zeros(T2, length(model.z))
        )
        for i in eachindex(angles)
    ]

    Q = zeros(T2, length(model.z))
    F = zeros(T2, length(model.z))
    χ = zeros(T2, length(model.z))
    S = zeros(T2, length(model.z))

    HeatingSolver(eos, opacities, angles, weights, bins, m, χ, S, Q, F, angle_solver)
end






#================================ Solving ====================================#

"""
    solve(h::HeatingSolver)

Solve the radiative heating equation with the given `solver`.
Frequency integration weights can be given via `weights`.
"""
function solve(h::HeatingSolver; weights=ones(h.bins))
    h.Q .= 0.0
    h.F .= 0.0

    for bin in 1:h.bins
        # Get the source function and opacity for this bin
        binoptics!(h, bin)

        # compute the radiation for all angles
        for (j, μ) in enumerate(h.angles)
            solve!(h.μ_solver[j], h.model.z, h.χ, h.S)
            h.Q .+= h.μ_solver[j].Q .* 
                    h.μ_solver[j].ω .*
                    h.χ .*
                    weights[bin]
            h.F .+= (h.μ_solver[j].Q .+ h.S).* 
                    h.μ_solver[j].ω .*
                    weights[bin]
        end
    end

    HeatingField(h.model, h.Q, h.F)
end

"""
    solve(h::AngleHeatingSolver, z, χ, S)

Solve the radiative heating equation in a given direction with given optics.
"""
function solve!(h::AngleHeatingSolver, z, χ, S)
    boundary!(h, z, χ, S)
    _solve!(
        h.μ,
        h.Q,
        h.dtau1, 
        h.dtau2, 
        h.ex0, 
        h.ex1, 
        h.ex2, 
        h.dsdtau1, 
        h.dsdtau2, 
        h.d1sdtau1, 
        h.d2sdtau2,
        h.χ, 
        h.S,
        h.z
    )
    boundary!(h)
end

"""
    _solve!(args...)

Qr solution from Dispatch (solver2, rt_mod.f90, Oct. 12, 2023)
"""
_solve!(
    μ,
    Q,
    dtau1, 
    dtau2, 
    ex0, 
    ex1, 
    ex2, 
    dsdtau1, 
    dsdtau2, 
    d1sdtau1, 
    d2sdtau2,
    χ, 
    S,
    z
) = begin
    aμ = abs(μ) 
    nz = length(z)
    @inbounds for i in 2:nz-1
        dtau1[i] = z[i] / aμ * 0.5 * (χ[i] + χ[i-1]) 
        dtau2[i] = z[i] / aμ * 0.5 * (χ[i] + χ[i+1])
        
        ex0[i] = exp(-dtau1[i])
        ex1[i] = 1 - ex0[i]
        ex2[i] = ex1[i] - dtau1[i] * ex0[i]
        dsdtau1[i] = (S[i] - S[i-1]) / dtau1[i]
        dsdtau2[i] = (S[i+1] - S[i]) / dtau2[i]
        d1sdtau1[i] = (dsdtau2[i]*dtau1[i] + dsdtau1[i]*dtau2[i])/(dtau1[i]+dtau2[i])
        d2sdtau2[i] = (dsdtau2[i] - dsdtau1[i])*2.0/(dtau1[i]+dtau2[i])
        
        Q[i] = Q[i-1]*ex0[i] - d1sdtau1[i] *ex1[i] + d2sdtau2[i] *ex2[i]
    end

    Q
end


"""
    binoptics!(solver, bin)

Compute source funtion and opacity for the given model in the given bin.
"""
function binoptics!(solver, bin)
    solver.S .= lookup(solver.eos, solver.opacities, :src, solver.model.lnρ, energyvariable(solver), bin)
    solver.χ .= lookup(solver.eos, solver.opacities, :κ,   solver.model.lnρ, energyvariable(solver), bin);
end

boundary!(h, z, χ, S) = begin
    if h.μ < 0
        h.z .= reverse(z)
        h.χ .= reverse(χ)
        h.S .= reverse(S)

        h.z[2:end] = abs.(diff(h.z)) / abs(h.μ)
        h.z[1] = h.z[2]
        
        dr = h.z[1]
        h.Q[1] = - h.S[1] * exp(-dr*h.χ[1])
        h.Q[end] = 0.0
    else
        h.z .= z
        h.χ .= χ
        h.S .= S

        h.z[2:end] = abs.(diff(h.z)) / abs(h.μ)
        h.z[1] = h.z[2]

        dr = h.z[end]
        h.Q[end] = - h.S[end] * exp(-dr*h.χ[end])
        h.Q[1] = 0.0
    end

    h
end

boundary!(h) = begin
    if h.μ < 0
        reverse!(h.Q)
    end

    h
end








#================================= Utilities =================================#

energyvariable(solver) = is_internal_energy(solver.eos) ? solver.model.lnEi : solver.model.lnT

heating(h::HeatingField) = h.values
flux(h::HeatingField) = h.flux
model(h::HeatingField) = h.model