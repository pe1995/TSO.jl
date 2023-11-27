#========================= Concrete types ====================================#

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

    Q_temp ::Vector{T}
    Q      ::Array{T, 2}
end

struct RadiationField{A<:AbstractModel, T<:AbstractFloat, N} <:AbstractRadiationField
    model::A
    values::Array{T, N}
end

struct InterpolatedRadiationField{R<:RadiationField, I} <:AbstractRadiationField
    radiation::R
    interpolation_function::Vector{I}
end






#========================== Constructors =====================================#

"""
Construct a binned radiative transfer solver based on a 1D model witht the kind:
    model[:, 1] -> z
    model[:, 2] -> ln(Ei)
    model[:, 3] -> ln(ρ)

    needs to be ordered from bottom to top. z increasing outwards.
"""
Solver(model::Array{T,2}; eos::AxedEoS, opacities::AbstractBinnedOpacities, angles=sort(labatto_4angles), weights=J_weights(labatto_4angles)) where {T} = begin
    opacities = opacities.opacities
    bins    = size(opacities.κ, 3)
    T2      = Float64 #eltype(opacities.κ)
    model   = Base.convert.(T2, deepcopy(model))

    # For lookup it is more convenient to save the log model
    model[:, 2] .= log.(model[:, 2])
    model[:, 3] .= log.(model[:, 3])

    if model[2, 3] > model[1, 3]
        # The second point has a higher density as the first one, so we are going down
        @warn "Model is increasing in density! Assuming flipped z axis..."
        model = reverse(model, dims=1)
    end

    if model[2, 1] < model[1, 1]
        @warn "z needs to be increasing outwards and start at the bottom. Flipping axes..."
        model[:, 1] .= -model[:, 1]
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
                        zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1)), zeros(T2, size(model, 1), bins))
end

Solver(model::AbstractModel, eos::AxedEoS; kwargs...) = begin
    model_array = zeros(eltype(model.z), length(model.z), 3)
    model_array[:, 1] = - deepcopy(model.z)
    model_array[:, 2] = exp.(getfield(model, energy_variable(eos)))
    model_array[:, 3] = exp.(model.lnρ)

    reverse!(model_array, dims=1)

    Solver(model_array; eos=eos, kwargs...)
end 

Solver(model::AbstractModel, eos::EoSTable; kwargs...) = Solver(model; eos=@axed(eos), kwargs...)

## For the inteprolation of the results
InterpolatedRadiationField(r::AbstractRadiationField, field::Symbol) = begin
    v = getfield(r.model, field)
    ip = [linear_interpolation(v |> reverse, r.values[:, i] .|> log |> reverse, extrapolation_bc=Line()) for i in axes(r.values, 2)]

    InterpolatedRadiationField(r, ip)
end

## Convenience
macro interpolated(r, field)
    r_l = esc(r)
    f_l = esc(field)

    :(InterpolatedRadiationField($(r_l), $(f_l)))
end





#=========================== Solving routines ================================#

function chat_transfer1d!(intensity, S, k, z, angles)
    n_points = length(z)

    source_function = similar(S)
    opacity = similar(k)
    depth   = similar(z)

    # Initialize the intensity array for each angle
    #intensity = zeros(n_points, length(angles))
    
    # Loop through each angle
    for (i, angle) in enumerate(angles)
        if angle < 0.0
            source_function .= reverse(S)
            opacity         .= reverse(k)
            depth           .= reverse(z)
            # Set the initial conditions for the intensity array for each angle
            intensity[1, i] = source_function[1]
        else
            source_function .= S
            opacity         .= k
            depth           .= z

            # Set the initial conditions for the intensity array for each angle
            intensity[1, i] = source_function[1]
        end

        # Loop through all the points, updating the intensity at each step
        for j in 2:n_points
            mean_opacity = 0.5 * (opacity[j-1] + opacity[j])
            intensity[j, i] = intensity[j-1, i] * exp(- mean_opacity * (depth[j] - depth[j-1]) / angle) + 
                               source_function[j-1] * (1.0 - exp(- mean_opacity * (depth[j] - depth[j-1]) / angle))
        end

        if angle < 0.0
            intensity[:, i] .= reverse(intensity[:, i])
        end

    end

    intensity
end

chat_transfer1d!(solver::BinnedTransferSolver) = chat_transfer1d!(solver.I, solver.S, solver.χ, solver.model[:, 1], solver.angles)

"""
Solve the 1D radiative transfer in a given direction μ=cos(θ).
"""
function transfer1d!(solver) 
    # Check in which direction we are going
    down = solver.μ[] .< 0.0

    @inbounds for i in axes(solver.model, 1)
        if down
            solver.z[i]      = solver.model[end+1-i, 1]
            solver.S_temp[i] = solver.S[end+1-i]
            solver.χ_temp[i] = solver.χ[end+1-i]
        else
            solver.z[i]      = solver.model[i, 1]
            solver.S_temp[i] = solver.S[i]
            solver.χ_temp[i] = solver.χ[i]
        end
    
        if i >= 2
            solver.dt[i-1] = 0.5 / solver.μ[] * (solver.z[i] - solver.z[i-1]) * ( solver.χ_temp[i] + solver.χ_temp[i-1] ) 
        end
    end

    S = solver.S_temp

    @inbounds for i in eachindex(solver.dt)
        solver.ex[i] = exp(-solver.dt[i]) 

        # Mask the small dt regions
        if solver.dt[i] > 1e-3
            solver.U0[i] = 1.0           - solver.ex[i]    
            solver.U1[i] = solver.dt[i]  - solver.U0[i]
        else
            solver.U1[i] = 0.5 * solver.dt[i]^2 - 1/6 * solver.dt[i]^2 * solver.dt[i]
            solver.U0[i] =  solver.dt[i] - solver.U1[i]
        end

        # Psi arrays
        solver.psiup[i] = solver.U0[i] - solver.U1[i] / solver.dt[i]
        solver.psi0[i]  = solver.U1[i] / solver.dt[i]

        # intensity at center point
        solver.I_temp[i] = S[i]

    end

    solver.I_temp[1] = first(S) *solver.ex[1]
    solver.Q_temp[1] = 0.0

    lI = length(solver.I_temp)

    # now solve
    for i in 2:lI
        solver.I_temp[i] = solver.I_temp[i-1]*solver.ex[i-1] + solver.psiup[i-1]*S[i-1] + solver.psi0[i-1]*S[i]
        #solver.Q_temp[i] = χ[i-1] *( solver.I_temp[i-1]*solver.ex[i-1] + solver.psiup[i-1]*S[i-1] + (solver.psi0[i-1]-1)*S[i])
        #@show i solver.I_temp[i] solver.psi0[i-1] solver.psiup[i-1] S[i]
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

        solver.Q[:, i] .= solver.Q_temp
        F[:, i] = mean_intensity ? J_from_I(solver) : flux_from_I(solver) 
    end

    F
end

function chat_binned_radiation(solver::B; mean_intensity=false) where {B<:BinnedTransferSolver}
    F = zeros(eltype(solver.I), size(solver.I, 1), solver.bins)

    #debug[] = true

    @inbounds for i in 1:solver.bins

        # set the source function and opacities
        optics!(solver, i)

        # solve the radiative transfer in this bin along all angles
        #for (j, angle) in enumerate(solver.angles)
        #    solver.μ[] = angle
        #    solver.I[:, j] .= transfer1d!(solver)
        #end

        chat_transfer1d!(solver)

        #if (!debug[])
        #    if debug_index[] == 0
        #        debug_index[] = i
        #    end
        #end
        solver.Q[:, i] .= solver.Q_temp
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





#=======================  Commputation of mean intensity =====================#

"""
Compute the mean intensity using the given Solver.
"""
Jν() = error("Please provide a solver.")

Jν(solver::BinnedTransferSolver) = sum(binned_radiation(solver, mean_intensity=true), dims=2)[:, 1]
Jν(solver::BinnedTransferSolver, weights) = begin
    f     = binned_radiation(solver, mean_intensity=true)
    mean_intensity(f, weights=weights)
end

chat_Jν(solver::BinnedTransferSolver) = sum(chat_binned_radiation(solver, mean_intensity=true), dims=2)[:, 1]
chat_Jν(solver::BinnedTransferSolver, weights) = begin
    f     = chat_binned_radiation(solver, mean_intensity=true)
    mean_intensity(f, weights=weights)
end

mean_intensity(I; weights=ones(size(I, 1))) = begin
    f_int = similar(I, size(I, 1))

    for i in axes(I, 1)
        f_int[i] = sum(weights .* I[i, :])
    end

    #solver.Q_temp .= I[1, :]
    f_int
end


"""
    mean_intensity(aos, opacity, model)

Compute the frequency dependent, angle integrated mean intensity from the given model.
Return a RadiationField object, that can be used for interpolation/extrapolation.
"""
mean_intensity(aos::AxedEoS, opacity::AbstractBinnedOpacities, model) = begin
    solver = Solver(model, aos, opacities=opacity)

    # Make a new model with the flipped axes saved
    args = []
    for f in fieldnames(typeof(model))
        v = if f == :lnρ
            solver.model[:, 3]
        elseif f == energy_variable(aos)
            solver.model[:, 2]
        elseif f == :z
            solver.model[:, 1]
        else
            getfield(model, f)
        end

        append!(args, [v])
    end

    m = typeof(model)(args...)

    RadiationField(m, binned_radiation(solver, mean_intensity=true))
end






#==================== Commputation of heating rate ===========================#

Qr(solver::BinnedTransferSolver) = begin
    I = binned_radiation(solver, mean_intensity=true)
    heating(solver, I)
end

Qr(solver::BinnedTransferSolver, weights) = begin
    I = binned_radiation(solver, mean_intensity=true)
    heating(solver, I; weights=weights)
end

chat_Qr(solver::BinnedTransferSolver, weights) = begin
    I = chat_binned_radiation(solver, mean_intensity=true)
    heating(solver, I; weights=weights)
end

chat_Qr(solver::BinnedTransferSolver) = begin
    I = chat_binned_radiation(solver, mean_intensity=true)
    heating(solver, I)
end

heating(solver, I; weights=ones(solver.bins)) = begin
    Q = zeros(eltype(I), size(I, 1))
    #S = similar(solver.S)
    for i in 1:solver.bins

        # set the source function and opacities
        optics!(solver, i)

        #S .= 0.0
        #for j in eachindex(solver.weights)
        #    S .+= solver.S .* solver.weights[j]
        #end
       
        #if i==1
        #    for j in axes(I, 1)
        #        @show i j I[j, i] solver.S[j] solver.χ[j]
        #    end
        #end

        #Q .+= solver.χ ./ exp.(solver.model[:, 3]) .* weights[i] .* (I[:, i] .- solver.S)
        Q .+= solver.χ .* weights[i] .* (I[:, i] .- solver.S)
    end

    Q .* 4 .*π 
end






#================ Source function integration for testing ====================#

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





#================== Opacity integration for testing ==========================#

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




#========================== Interpolation API ================================#

evaluate_radiation(r::InterpolatedRadiationField, bin_index, values...) = r.interpolation_function[bin_index](values...)
broadcastable(r::InterpolatedRadiationField) = Ref(r)






#======================== General methods ====================================#

effective_temperature(F) = (F/σ_S)^(1/4)
@inline angle_index(solver) = findfirst(solver.angles .≈ solver.μ[])
J_weights(angles) = ω_midpoint(sort(angles), -1, 1) ./ 2





#====================== General global variables =============================#
const debug = Ref(true)
const debug_index = Ref(0)

const labatto_4angles  = [0.9739065285, 0.8650633666, 0.6794095682, 0.4333953941, -0.9739065285, -0.8650633666, -0.6794095682, -0.4333953941]
const labatto_4weights = [0.0666713443, 0.1494513491, 0.2190863625, 0.2692667193,  0.0666713443,  0.1494513491,  0.2190863625,  0.2692667193]
const labatto_7angles  = [-0.9951872199, -0.9604912684, -0.8884592328, -0.7818314825, -0.6423493394, -0.4646021794, -0.2492867453, 
                           0.2492867453,  0.4646021794,  0.6423493394,  0.7818314825,  0.8884592328,  0.9604912684,  0.9951872199]