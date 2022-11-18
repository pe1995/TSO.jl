abstract type RadiativeTransferSolver end

struct BinnedTransferSolver{T<:AbstractFloat} <:RadiativeTransferSolver
    eos       ::EoSTable
    opacities ::OpacityTable
    angles    ::Vector{T}
    weights   ::Vector{T}
    bins      ::Int
    model     ::Array{T, 2}
    χ         ::Vector{T}
    S         ::Vector{T}
    I         ::Array{T, 2}
end


## Constructors

"""
Construct a binned radiative transfer solver based on a 1D model witht the kind:
    model[:, 1] -> z
    model[:, 2] -> ln(Ei)
    model[:, 3] -> ln(ρ)
"""
Solver(model::Array{T,2}; eos, opacities, angles=[1.0, 0.76505530, 0.28523153], weights=[0.13333334, 0.14477713*4, 0.07913155*4]) where {T} = begin
    bins    = size(opacities.κ, 3)
    T2      = eltype(opacities.κ)
    model   = Base.convert.(T2, model)

    # For lookup it is more convenient to save the log model
    model[:, 2] .= log.(model[:, 2])
    model[:, 3] .= log.(model[:, 3])

    angles  = Base.convert.(T2, angles)
    weights = Base.convert.(T2, weights)
    BinnedTransferSolver(eos, opacities, angles, weights, bins, model, zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1), length(angles)))
end


## Solving routines

"""Solve the 1D radiative transfer in a given direction μ=cos(θ)."""
function transfer1d(solver, μ)
    z = view(solver.model, :, 1)

    dt   = 0.5 .* μ .* abs.(z[2:end] .- z[1:end-1]) .* ( solver.χ[2:end] + solver.χ[1:end-1] ) 
    dt2  = dt .^2
    ex   = exp.(-dt) 

    # Compute the U functions
    U0 = 1.0 .- ex
    U1 = dt  .- U0

    psiup = U0 .- U1 ./ dt
    psi0  = U1 ./ dt

    # intensity at center point
    I    = similar(solver.S)
    I[1] = first(solver.S)
    lI   = length(I)

    for i in 2:lI
        I[i] = I[i-1]*ex[i-1] + psiup[i-1]*solver.S[i-1] + psi0[i-1]*solver.S[i]
    end

    I
end

"""Integrate the angles of I."""
flux_from_I(solver) = begin
    f = similar(solver.I, size(solver.I, 1))
    for i in axes(solver.I, 1)
        f[i] = sum(solver.weights .* solver.angles .* solver.I[i, :] *2π)
    end

    f
end

J_from_I(solver) = begin
    f = similar(solver.I, size(solver.I, 1))
    for i in axes(solver.I, 1)
        f[i] = 1/(4π) *sum(solver.weights .* solver.I[i, :] *2π)
    end

    f
end

"""Compute source funtion and opacity for the given model in the given bin."""
function optics!(solver, bin)
    solver.S .= exp.(lookup(solver.eos, solver.opacities, :src, view(solver.model, :, 3), view(solver.model, :, 2), bin))
    solver.χ .= exp.(lookup(solver.eos, solver.opacities, :κ,   view(solver.model, :, 3), view(solver.model, :, 2), bin));
end

"""Compute the flux from the given solver from Binned radiative transfer."""
function binned_radiation(solver::BinnedTransferSolver; mean_intensity=false)
    F = zeros(eltype(solver.I), size(solver.I, 1), solver.bins)

    for i in 1:solver.bins

        # set the source function and opacities
        optics!(solver, i)

        # solve the radiative transfer in this bin along all angles
        for (j, angle) in enumerate(solver.angles)
            solver.I[:, j] .= transfer1d(solver, angle)
        end

        F[:, i] = mean_intensity ? J_from_I(solver) : flux_from_I(solver) 
    end

    F
end

flux(solver::BinnedTransferSolver) = sum(binned_radiation(solver), dims=2)
flux(solver::BinnedTransferSolver, weights) = begin
    f     = binned_radiation(solver)
    f_int = similar(f, size(f, 1))

    for i in axes(f, 1)
        f_int[i] = sum(weights .* f[i, :])
    end

    f_int
end



## Commputation of mean intensity

Jν(solver::BinnedTransferSolver) = sum(binned_radiation(solver, mean_intensity=true), dims=2)
Jν(solver::BinnedTransferSolver, weights) = begin
    f     = binned_radiation(solver, mean_intensity=true)
    f_int = similar(f, size(f, 1))

    for i in axes(f, 1)
        f_int[i] = sum(weights .* f[i, :])
    end

    f_int
end



## General methods

effective_temperature(F) = (F/σ_S)^(1/4)

