abstract type RadiativeTransferSolver end

struct BinnedTransferSolver{T<:AbstractFloat} <:RadiativeTransferSolver
    eos       ::EoSTable
    opacities ::OpacityTable
    angles    ::T
    weights   ::T
    bins      ::T
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
BinnedTransferSolver(model::Array{T,2}; eos, opacities, angles=[1.0, 0.76505530, 0.28523153], weights=[0.13333334, 0.14477713*4, 0.07913155*4]) where {T} = begin
    bins    = size(opacities.κ, 3)
    T2      = eltype(opacities.κ)
    model   = Base.convert.(T2, model)
    angles  = Base.convert.(T2, angles)
    weights = Base.convert.(T2, weights)
    BinnedTransferSolver(eos, opacities, angles, weights, bins, model, zeros(size(model, 1)), zeros(size(model, 1)), zeros(size(model, 1)))
end


## Solving routines

"""Solve the 1D radiative transfer in a given direction μ=cos(θ)."""
function transfer1d(solver, μ)
    dt   = 0.5 .* μ .* (solver.z[2:end] .- solver.z[1:end-1]) .* ( solver.χ[2:end] + solver.χ[1:end-1] ) 
    dt2  = dt .^2
    ex   = exp.(-dt) 

    # Compute the U functions
    U0 = 1.0 .- ex
    U1 = dt  .- U0

    psiup = U0 .- U1 ./ dt
    psi0  = U1 ./ dt

    # intensity at center point
    
    I    = similar(solver.z)
    I[1] = first(solver.S)

    for i in eachindex(view(I, 2:length(I)))
        I[i] = solver.I[i-1, iμ] *ex[i-1] + psiup[i-1] * solver.S[i-1] + psi0[i-1] * solver.S[i]
    end

    I
end

"""Integrate the angles of I."""
flux_from_I(solver) = sum(solver.weights .* solver.I, dims=2)

"""Compute source funtion and opacity for the given model in the given bin."""
function optics!(solver, bin)
    solver.S .= lookup(solver.eos, solver.opacites, :src, view(solver.model, :, 3), view(solver.model, :, 2), bin)
    solver.χ .= lookup(solver.eos, solver.opacites, :κ,   view(solver.model, :, 3), view(solver.model, :, 2), bin);
end

"""Compute the flux from the given solver from Binned radiative transfer."""
function flux(solver::BinnedTransferSolver)
    F = zeros(length(solver.I), length(solver.bins))

    for i in enumerate(solver.bins)

        # set the source function and opacities
        optics!(solver, i)

        # solve the radiative transfer in this bin along all angles
        for (j,angle) in enumerate(solver.angles)
            solver.I[:, j] .= transfer1d(solver, angle)
        end

        F[:, i] = last(flux_from_I(solver))
    end

    sum(F, dims=2)
end

