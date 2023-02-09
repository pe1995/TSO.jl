## Concrete types

struct BinnedTransferSolver{T<:AbstractFloat, A<:AxedEoS, O<:OpacityTable} <:RadiativeTransferSolver
    eos       ::A
    opacities ::O
    angles    ::Vector{T}
    weights   ::Vector{T}
    bins      ::Int
    model     ::Array{T, 2}
    χ         ::Vector{T}
    S         ::Vector{T}
    I         ::Array{T, 2}
    
    ## Process variables
    μ      ::Ref{T}
    z      ::Vector{T}     
    dt     ::Vector{T}
    ex     ::Vector{T}
    U0     ::Vector{T}
    U1     ::Vector{T}
    mask   ::BitVector
    psiup  ::Vector{T}
    psi0   ::Vector{T}
    I_temp ::Vector{T}
    S_temp ::Vector{T}
    χ_temp ::Vector{T}
end




## Constructors

"""
Construct a binned radiative transfer solver based on a 1D model witht the kind:
    model[:, 1] -> z
    model[:, 2] -> ln(Ei)
    model[:, 3] -> ln(ρ)

    needs to be ordered from bottom to top. z increasing outwards.
"""
Solver(model::Array{T,2}; eos::AxedEoS, opacities::BinnedOpacities, angles=[1.0, 0.76505530, 0.28523153, -1, -0.76505530, -0.28523153], weights=[0.13333334/2, 0.14477713*2, 0.07913155*2, 0.13333334/2, 0.14477713*2, 0.07913155*2]) where {T} = begin
    opacities = opacities.opacities
    bins    = size(opacities.κ, 3)
    T2      = eltype(opacities.κ)
    model   = Base.convert.(T2, deepcopy(model))

    # For lookup it is more convenient to save the log model
    model[:, 2] .= log.(model[:, 2])
    model[:, 3] .= log.(model[:, 3])

    if model[2, 1] < model[1,1]
        @warn "z needs to be increasing outwards and start at the bottom. Flipping axes..."
        model = reverse(model, dims=1)
    end

    mask    = sortperm(angles, rev=true)
    angles  = Base.convert.(T2, angles[mask])
    weights = Base.convert.(T2, weights[mask])

    BinnedTransferSolver(eos, opacities, angles, weights, bins, model, 
                        zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1), length(angles)), 
                        Ref{T2}(T2(0.0)),
                        zeros(T2, size(model, 1)), 
                        zeros(T2, size(model, 1)-1), zeros(T2, size(model, 1)-1), zeros(T2, size(model, 1)-1), zeros(T2, size(model, 1)-1),
                        falses(size(model, 1)-1),    
                        zeros(T2, size(model, 1)-1), zeros(T2, size(model, 1)-1), 
                        zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1)))
end

Solver(model::Model1D, eos::AxedEoS; kwargs...) = begin
    model_array = zeros(eltype(model.z), length(model.z), 3)
    model_array[:, 1] = deepcopy(model.z)
    model_array[:, 2] = exp.(getfield(model, energy_variable(eos)))
    model_array[:, 3] = exp.(model.lnρ)

    Solver(model_array; eos=eos, kwargs...)
end 

Solver(model::Model1D, eos::EoSTable; kwargs...) = Solver(model; eos=@axed(eos), kwargs...)



## Solving routines

"""
Solve the 1D radiative transfer in a given direction μ=cos(θ).
"""
function transfer1d!(solver) 

    # Check in which direction we are going
    down = solver.μ[] .< 0.0

    solver.z      .= solver.model[:, 1]
    solver.S_temp .= solver.S
    solver.χ_temp .= solver.χ

    if down
        solver.z      .= reverse(solver.z) 
        solver.S_temp .= reverse(solver.S)
        solver.χ_temp .= reverse(solver.χ)
    end

    S = solver.S_temp
    χ = solver.χ_temp
    z = solver.z

    solver.dt .= 0.5 .* solver.μ[] .* (z[2:end] .- z[1:end-1]) .* ( χ[2:end] + χ[1:end-1] ) 

    #@assert all(solver.dt .> 0.0)

    @inbounds for i in eachindex(solver.dt)
        solver.ex[i] = exp(-solver.dt[i]) 

        # Compute the U functions
        solver.U0[i] = 1.0           - solver.ex[i]
        solver.U1[i] = solver.dt[i]  - solver.U0[i]

        # Mask the small dt regions
        if solver.dt[i] < 1e-3
            solver.U1[i] = 0.5 * solver.dt[i]^2 - 1/6 * solver.dt[i]^2 * solver.dt[i]
            solver.U0[i] =  solver.dt[i] - solver.U1[i]
        end

        # Psi arrays
        solver.psiup[i] = solver.U0[i] - solver.U1[i] / solver.dt[i]
        solver.psi0[i]  = solver.U1[i] / solver.dt[i]

        # intensity at center point
        solver.I_temp[i] = S[i]
    end

    solver.I_temp[1] = if !down
        first(S) 
    else
        solver.psiup[1]*S[1] + solver.psi0[1]*S[1]
        first(S) 
    end 

    lI = length(solver.I_temp)

    # now solve
    for i in 2:lI
        solver.I_temp[i] = solver.I_temp[i-1]*solver.ex[i-1] + solver.psiup[i-1]*S[i-1] + solver.psi0[i-1]*S[i]
    end

    down ? reverse(solver.I_temp) : solver.I_temp
end

"""
Integrate the angles of I.
"""
flux_from_I(solver) = begin
    f = similar(solver.I, size(solver.I, 1))
    for i in axes(solver.I, 1)
        f[i] = sum(solver.weights .* solver.angles .* solver.I[i, :])
    end

    f #./ sum(solver.weights)
end

J_from_I(solver) = begin
    f = similar(solver.I, size(solver.I, 1))
    for i in axes(solver.I, 1)
        f[i] = sum(solver.weights .* solver.I[i, :])
    end

    f #./ sum(solver.weights)
end

"""
Compute source funtion and opacity for the given model in the given bin.
"""
function optics!(solver, bin)
    solver.S .= lookup(solver.eos, solver.opacities, :src, view(solver.model, :, 3), view(solver.model, :, 2), bin)
    solver.χ .= lookup(solver.eos, solver.opacities, :κ,   view(solver.model, :, 3), view(solver.model, :, 2), bin);
end

"""
Compute the radiation field from the given solver from Binned radiative transfer.
"""
function binned_radiation(solver::B; mean_intensity=false) where {B<:BinnedTransferSolver}
    F = zeros(eltype(solver.I), size(solver.I, 1), solver.bins)

    #debug[] = true

    @inbounds for i in 1:solver.bins

        # set the source function and opacities
        optics!(solver, i)

        # solve the radiative transfer in this bin along all angles
        for (j, angle) in enumerate(solver.angles)
            solver.μ[] = angle
            solver.I[:, j] .= transfer1d!(solver)
        end

        #if (!debug[])
        #    if debug_index[] == 0
        #        debug_index[] = i
        #    end
        #end

        F[:, i] = mean_intensity ? J_from_I(solver) : flux_from_I(solver) 
    end

    F
end

"""
Compute the Flux using the given Solver.
"""
flux() = error("Please provide a solver.")

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

"""
Compute the mean intensity using the given Solver.
"""
Jν() = error("Please provide a solver.")

Jν(solver::BinnedTransferSolver) = sum(binned_radiation(solver, mean_intensity=true), dims=2)
Jν(solver::BinnedTransferSolver, weights) = begin
    f     = binned_radiation(solver, mean_intensity=true)
    mean_intensity(f, weights=weights)
end

mean_intensity(I; weights=ones(size(I, 1))) = begin
    f_int = similar(I, size(I, 1))

    for i in axes(I, 1)
        f_int[i] = sum(weights .* I[i, :])
    end

    f_int
end





## Commputation of heating rate

Qr(solver::BinnedTransferSolver) = begin
    I = binned_radiation(solver, mean_intensity=true)
    heating(solver, I)
end

Qr(solver::BinnedTransferSolver, weights) = begin
    I = binned_radiation(solver, mean_intensity=true)
    heating(solver, I; weights=weights)
end

heating(solver, I; weights=ones(solver.bins)) = begin
    Q = zeros(eltype(I), size(I, 1))
    S = similar(solver.S)
    for i in 1:solver.bins

        # set the source function and opacities
        optics!(solver, i)

        #S .= 0.0
        #for j in eachindex(solver.weights)
        #    S .+= solver.S .* solver.weights[j]
        #end
        Q .+= solver.χ .* weights[i] .* (I[:, i] .- S)
    end

    Q .* 4 .*π 
end





## Source function integration for testing

S(solver::BinnedTransferSolver) = begin
    S_b = zeros(eltype(solver.S), length(solver.S))
    source_function(solver, S_b)
end

S(solver::BinnedTransferSolver, weights) = begin
    S_b = zeros(eltype(solver.S), length(solver.S))
    source_function(solver, S_b, weights=weights)
end

source_function(solver::BinnedTransferSolver, S_b; weights=ones(solver.bins)) = begin
    for i in 1:solver.bins
        # set the source function and opacities in this bin
        optics!(solver, i)

        for j in eachindex(solver.weights)
            S_b .+= weights[i] * solver.weights[j] .* solver.S
        end
    end

    S_b
end





## Opacity integration for testing

κ(solver::BinnedTransferSolver) = begin
    S_b = zeros(eltype(solver.S), length(solver.S))
    opacity(solver, S_b)
end

opacity(solver::BinnedTransferSolver, S_b) = begin
    for i in 1:solver.bins
        # set the source function and opacities in this bin
        optics!(solver, i)

        S_b .+= solver.χ
    end

    S_b / solver.bins
end





## General methods

effective_temperature(F) = (F/σ_S)^(1/4)
@inline angle_index(solver) = findfirst(solver.angles .≈ solver.μ[])





## General global variables
const debug = Ref(true)
const debug_index = Ref(0)