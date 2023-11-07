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
Average3D(eos, path; logg=log10(2.75e4))      = begin
    model = Average3D(path)
    e     = lnT_to_lnEi(eos, model.lnρ, model.lnT)
    Model1D(model.z, model.lnρ, model.lnT, e, logg)
end
Average3D(path; logg=log10(2.75e4)) = begin
    solar_model        = reverse(readdlm(path, skipstart=0), dims=1)
    solar_model[:, 1] .= -solar_model[:, 1]
    solar_model[:, 3] .= exp.(solar_model[:, 3])
    z, lnρ, lnT        = solar_model[:, 1], log.(solar_model[:, 3]), log.(solar_model[:, 2]) 
    Model1D(z, lnρ, lnT, zeros(length(z)), logg)
end

macro optical(model, eos, opacities)
    m = esc(model)
    e = esc(eos)
    o = esc(opacities)
    quote
        OpticalModel1D(Base.convert(typeof(getfield($m, :z)), rosseland_optical_depth($e, $o, $m)), (getfield($m, p) for p in fieldnames(typeof($m)))...)
    end
end

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
        Model1D(input_arrays[2:end]..., logg)
    else
        OpticalModel1D(input_arrays..., logg)
    end
end



#======================================================== Model conversions ==#

"""
Convert z,T,ρ to z,E,ρ model based on EoS.
"""
lnT_to_lnEi(eos::AxedEoS, lnρ, lnT) = lookup(eos, :lnEi, lnρ, lnT)
lnT_to_lnEi(eos::RegularEoSTable, lnρ, lnT) = lookup(AxedEoS(eos), :lnEi, lnρ, lnT)



#================================================================ Utilities ==#

upsample(model::AbstractModel, N=500) = begin
	oldz = model.z
	z = Base.convert.(eltype(oldz), 
		range(minimum(oldz), maximum(oldz), length=N)) |> collect

	is_field(model, f) = length(getfield(model, f)) == length(oldz)
	fields = [f for f in fieldnames(typeof(model)) if is_field(model, f)]

	
	results = Dict(f=>similar(model.lnT, N) for f in fields)
	
	for (i, f) in enumerate(fields)
		v = getfield(model, f)
		results[f] .= TSO.linear_interpolation(oldz, v).(z)
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
It used the density to determine where the bottom is. It will fo from bottom to top and 
with the lowest z value at the bottom.
"""
function flip!(model::AbstractModel)
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

    model
end

function flip(model::AbstractModel)
    m = deepcopy(model)
    flip!(m)

    m
end

function pick_point(model, i)
	args = Dict()
	for f in fieldnames(typeof(model))
		v = getfield(model, f)
		if typeof(v) <: AbstractArray
			args[f] = [v[i]]
		else
			args[f] = v
		end
	end

	Model1D(;args...)
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

#=============================================================================#
