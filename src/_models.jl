#================================================= Model Atmosphere structs ==#

struct Model1D{T} <:AbstractModel
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
    logg ::T
end


struct OpticalModel1D{T} <:AbstractModel
    τ    ::Vector{T}
    z    ::Vector{T}
    lnρ  ::Vector{T}
    lnT  ::Vector{T}
    lnEi ::Vector{T}
    logg ::T
end



#============================================================= Constructors ==#

Model1D(eos, z, lnρ, lnT; logg=log10(2.75e4)) = Model1D(z, lnρ, lnT, lnT_to_lnEi(eos, lnρ, lnT), logg)

"""
    Average3D(eos, path; logg=log10(2.75e4))

Construct a Model1D from 3 or 4 column delimited file. Add `lnEi` from given EoS.
"""
Average3D(eos, path; logg=log10(2.75e4)) = begin
    model = Average3D(path, logg=logg)
    e = lnT_to_lnEi(eos, model.lnρ, model.lnT)

    paras = Dict(f=>getfield(model, f) for f in fieldnames(typeof(model)))
    paras[:lnEi] .= e

    Model1D(; paras...)
end

"""
    Average3D(path; logg=log10(2.75e4))

Construct a Model1D from 3 or 4 column delimited file.
"""
Average3D(path; logg=log10(2.75e4)) = begin
    solar_model        = reverse(readdlm(path, skipstart=0), dims=1)
    solar_model[:, 1] .= -solar_model[:, 1]
    solar_model[:, 3] .= exp.(solar_model[:, 3])
    z, lnρ, lnT        = solar_model[:, 1], log.(solar_model[:, 3]), log.(solar_model[:, 2]) 

    τ = if size(solar_model, 2) == 4
        exp.(solar_model[:, 4])
    else
        nothing 
    end

    Model1D(z=z, lnρ=lnρ, lnT=lnT, lnEi=zeros(length(z)), logg=logg, τ=τ)
end

macro optical(model, eos, opacities)
    m = esc(model)
    e = esc(eos)
    o = esc(opacities)
    quote
        OpticalModel1D(Base.convert(typeof(getfield($m, :z)), rosseland_optical_depth($e, $o, $m)), (getfield($m, p) for p in fieldnames(typeof($m)))...)
    end
end

macro optical(model, eos)
    m = esc(model)
    e = esc(eos)
    quote
        OpticalModel1D(Base.convert(typeof(getfield($m, :z)), rosseland_optical_depth($e, $m)), (getfield($m, p) for p in fieldnames(typeof($m)) if p!=:τ)...)
    end
end


"""
    Model1D(; kwargs...)

Construct 1D model from keyword arrays. If `τ` given OpticalModel1D is constructed.

# Examples

```julia
julia> Model1D(τ=nothing, z=nothing, lnρ=nothing, lnT=nothing, lnEi=nothing, logg=log10(2.75e4))
TSO.OpticalModel1D{Float64}(...)
```
"""
Model1D(; τ   =nothing,
          z   =nothing,
          lnρ =nothing,
          lnT =nothing,
          lnEi=nothing,
          logg=log10(2.75e4)) = begin
    
    model_arrays = [τ, z, lnρ, lnT, lnEi]
    names        = [:τ, :z, :lnρ, :lnT, :lnEi]

    @assert any(.!isnothing.(model_arrays))

    i_ref = findfirst(.!isnothing.(model_arrays))
    model_arrays = if length(size(model_arrays[i_ref]))==0
        model_arrays = [[i] for i in model_arrays]
    else
        model_arrays
    end

    s_ref = size(model_arrays[i_ref])
    t_ref = eltype(model_arrays[i_ref])

    input_arrays = [zeros(t_ref, s_ref...) for _ in names]
    for (i,para) in enumerate(names)
        if isnothing(model_arrays[i])
            continue
        elseif isnothing(first(model_arrays[i]))
            continue
        else
            input_arrays[i] = Base.convert.(t_ref, model_arrays[i])
        end
    end

    if isnothing(τ)
        Model1D(input_arrays[2:end]..., Base.convert(t_ref, logg))
    else
        OpticalModel1D(input_arrays..., Base.convert(t_ref, logg))
    end
end



#======================================================== Model conversions ==#

"""
    lnT_to_lnEi(eos::AxedEoS, lnρ, lnT)

Convert z,T,ρ to z,E,ρ model based on EoS.
"""
lnT_to_lnEi(eos::AxedEoS, lnρ, lnT) = lookup(eos, :lnEi, lnρ, lnT)
lnT_to_lnEi(eos::RegularEoSTable, lnρ, lnT) = lookup(AxedEoS(eos), :lnEi, lnρ, lnT)


#=================================================================== optics ==#

"""
    optical_surface(model)

Return the z coordinate of the optical surface
"""
optical_surface(model::OpticalModel1D) = optical_surface(model.τ, model.z)

"""
    optical_surface!(model)

Move the 0 point to the optical surface.
"""
optical_surface!(model::OpticalModel1D) = model.z .= model.z .- optical_surface(model)

optical_surface(tau, z) = begin
    mask = sortperm(log10.(tau))
    linear_interpolation(
        Interpolations.deduplicate_knots!(log10.(tau[mask]), move_knots=true), 
        z[mask],
        extrapolation_bc=Line()
    )(0.0)
end

#================================================================ Utilities ==#

upsample(model::AbstractModel, N=500) = begin
	oldz = deepcopy(model.z)
	z = Base.convert.(eltype(oldz), 
		range(minimum(oldz), maximum(oldz), length=N)) |> collect

	is_field(model, f) = length(getfield(model, f)) == length(oldz)
	fields = [f for f in fieldnames(typeof(model)) if is_field(model, f)]

	
	results = Dict(f=>similar(model.lnT, N) for f in fields)
	
	for (i, f) in enumerate(fields)
		v = getfield(model, f)
		results[f] .= TSO.linear_interpolation(Interpolations.deduplicate_knots!(oldz, move_knots=true), v).(z)
	end

	args = [!is_field(model, f) ? getfield(model, f) : results[f] 
				for f in fieldnames(typeof(model))]

	typeof(model)(args...)
end

function cut(m::AbstractModel; kwargs...)
	masks = trues(length(m.z))
	for (para, lims) in kwargs
		v = getfield(m, para)
		mask = first(lims) .< v .< last(lims)
		masks .= masks .& mask
	end

	dat = []
	for f in fieldnames(typeof(m))
		v = getfield(m, f)
		if typeof(v) <: AbstractArray
			append!(dat, [v[masks]])
		else
			append!(dat, [v])
		end
	end

	typeof(m)(dat...)
end

function interpolate_to(m::AbstractModel; in_log=false, kwargs...)
    @assert length(keys(kwargs)) == 1
    
    para, new_scale = first(keys(kwargs)), first(values(kwargs))
    old_scale = in_log ? log10.(getfield(m, para)) : getfield(m, para)
    old_scale = deepcopy(old_scale)
    mask = sortperm(old_scale)

    d = []
    for f in fieldnames(typeof(m))
        if f == para
            in_log ? append!(d, [exp10.(new_scale)]) : append!(d, [new_scale])
        else
            v = getfield(m, f)
            if typeof(v) <: AbstractArray
                r = linear_interpolation(
                    Interpolations.deduplicate_knots!(old_scale[mask], move_knots=true),
                    v[mask], 
                    extrapolation_bc=Line()
                ).(new_scale)

                # the upper layers are almost isothermal
                if f==:lnT
                    mint = minimum(v[mask])
                    r[r.<mint] .= mint
                end

                append!(d, [r])
            else
                append!(d, [v])
            end
        end
    end

    flip!(typeof(m)(d...))
end


"""
    flip!(model)

Flip the box object and possibly reverse the sign of the height scale (only geometrical).
It used the density to determine where the bottom is. It will be from bottom to top and 
with the lowest z value at the bottom. `depth` will be the other way around.
"""
function flip!(model::AbstractModel; depth=false)
    d = exp.(model.lnρ)
    is_upside_down = first(d) < last(d)

    if is_upside_down
        for f in fieldnames(typeof(model))
            v = getfield(model, f)
            if typeof(v) <: AbstractArray
                reverse!(v)
            end
        end
    end

    # Now it is from bottom to top, so the first value in z should be the smalles value
    if model.z[1] > model.z[end]
        model.z .*= -1
    end

    # if depth, make the opposite now
    if depth
        # it is bottom to top
        for f in fieldnames(typeof(model))
            v = getfield(model, f)
            if typeof(v) <: AbstractArray
                reverse!(v)
            end
        end

        if model.z[1] > model.z[end]
            model.z .*= -1
        end
    end

    model
end

"""
    flip(model)

Flip the box object and possibly reverse the sign of the height scale (only geometrical).
It used the density to determine where the bottom is. It will be from bottom to top and 
with the lowest z value at the bottom. `depth` will be the other way around.
"""
function flip(model::AbstractModel; kwargs...)
    m = deepcopy(model)
    flip!(m; kwargs...)

    m
end

_pick_point(model, pick_index) = begin
    args = Dict()
	for f in fieldnames(typeof(model))
		v = getfield(model, f)
		if typeof(v) <: AbstractArray
			args[f] = [pick_index(v)]
		else
			args[f] = v
		end
	end

	Model1D(; args...)
end

function pick_point(model, i)
    pick_index = (typeof(i)<:Function) ? i : a -> a[i]
	_pick_point(model, pick_index)
end

function convert_model(model::AbstractModel, tp)
    flds = []
    for f in fieldnames(typeof(model))
        v = getfield(model, f)
        if typeof(v) <: AbstractArray
            append!(flds, [Base.convert.(tp, v)])
        else
            append!(flds, [Base.convert(tp, v)])
        end
    end

    typeof(model)(flds...)
end

Base.length(m::AbstractModel) = length(m.z)
Base.size(m::AbstractModel) = size(m.z)


#======================================== Make models monotonic in T and rho =#

"""
    monotonic(model; spline=true)

Make model monotonically increasing in rho and T. Sort it to be on a depth scale first, 
followed by a reverse minima accumulation and spline interpolation (SteffenMonotonicInterpolation), if wanted.
"""
monotonic(model; spline=true) = begin
    model_depth = flip(deepcopy(model), depth=true)
    x = copy(model_depth.z)
	y = similar(model_depth.lnρ)
    y .= reverse(model_depth.lnρ)
    y .= reverse(accumulate(min, y))

    model_depth.lnρ .= if spline
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

    y .= reverse(model_depth.lnT)
    y .= reverse(accumulate(min, y))
    model_depth.lnT .= if spline
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

    model_depth
end


#=============================================================================#

"""
    save_text_m3d(model_in::AbstractModel, args...; kwargs...)

Save the model in a format readable by M3D.
"""
save_text_m3d(model_in::AbstractModel, args...; kwargs...) = begin
    model = flip(model_in, depth=true)
    save_text_m3d(model.z, exp.(model.lnT), exp.(model.lnρ), args...; kwargs...)
end

"""
    save_text_m3d(z, T, ρ, f_new; header=nothing, vmic=zeros(length(z)))

Save the model in a format readable by M3D.
"""
function save_text_m3d(z, T, ρ, f_new; header=nothing, vmic=zeros(length(z)))
    open(f_new, "w") do f
        h = isnothing(header) ? "TSO.Model1D" : header
		write(f, h*"\n")
		write(f, "$(length(z))\n")
		
		for i in eachindex(z)
			line = @sprintf "%.6E    %.1f    %.4E    %.4E    %.1f\n" z[i] T[i] 0.0 ρ[i] vmic[i]
			write(f, line)
		end
	end
end
