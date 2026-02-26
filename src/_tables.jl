"""
Table containing binned Opacities + Source functions on a regular grid.
"""
struct RegularOpacityTable{T<:AbstractFloat, 
                            NOp, NRoss, NSrc} <: AbstractRegularTable
    κ             ::Array{T, NOp}
    κ_ross        ::Array{T, NRoss}
    src           ::Array{T, NSrc}
    λ             ::Array{T, 1}
    optical_depth ::Bool
end

"""
Table containing EoS properties on a regular grid.
"""
struct RegularEoSTable{F<:AbstractFloat, 
                        NRho, NT, NEi, 
                        NPg, NRoss, NNe} <: AbstractRegularTable
    lnRho  :: Array{F, NRho}
    lnT    :: Array{F, NT}
    lnEi   :: Array{F, NEi}
    lnPg   :: Array{F, NPg}
    lnRoss :: Array{F, NRoss}
    lnNe   :: Array{F, NNe}
end

"""
Table containing binned Opacities + Source functions on an irregular grid.
"""
struct IrregularOpacityTable <: AbstractIrregularTable
end

"""
Table containing EoS properties on an irregular grid.
"""
struct IrregularEoSTable <: AbstractIrregularTable
end

"""
Axis info type.
"""
struct EoSAxis{F<:AbstractFloat, Naxis, I<:Integer, R<:AbstractRange}
    values     ::Array{F, Naxis}
    dimension  ::I
    length     ::I
    name       ::Symbol
    unirange   ::R
    is_uniform ::Bool
end

"""
EoS Table with clearly defined Axes. Convenient for faster lookup if axes are uniform.
"""
struct AxedEoS{E<:RegularEoSTable, DA<:EoSAxis, EA<:EoSAxis}
    eos          ::E
    energy_axes  ::EA
    density_axes ::DA
end




#= Lookup functionality =# 

struct UniformLookup{F} <:LookupFunction
    f::F
end

struct GriddedLookup{F} <:LookupFunction
    f::F
end

struct ScatteredLookup{F} <:LookupFunction
    f::F
end

struct BisectionLookup{F} <:LookupFunction
    f::F
end

struct ZeroLookup{F} <:LookupFunction
    f::F
end




#= Aliases =#

OpacityTable = Union{RegularOpacityTable, IrregularOpacityTable} 
EoSTable     = Union{RegularEoSTable, IrregularEoSTable} 
SqOpacity    = RegularOpacityTable
SqEoS        = RegularEoSTable




#= Convenience Constructor functions =#

Opacity(args...; regular=true, kwargs...) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...;     regular=true, kwargs...) = regular ? RegularEoSTable(args...;     kwargs...) : IrregularEoSTable(args...;     kwargs...)

RegularOpacityTable(κ::Array, κ_ross::Array, src::Array; optical_depth=false) = RegularOpacityTable(κ, κ_ross, src, optical_depth) 
#RegularOpacityTable(; κ, src, κ_ross=zeros_as(κ, 1:3), κ_ross=zeros_as(κ, 1:3))

AxedEoS(eos) = AxedEoS(eos, EnergyAxis(eos), DensityAxis(eos))

"""
    @axed(eos)

Add axes information to the given EoS table.
"""
macro axed(eos)
    eos_l = esc(eos)
    :(AxedEoS($(eos_l)))
end

opacity(opa::RegularOpacityTable) = opa
wavelength(opa::RegularOpacityTable) = opa.λ
table(eos::AxedEoS) = table(eos.eos)
table(eos::RegularEoSTable) = eos




#= General functions =#

limits(aos::AxedEoS, args...; kwargs...) = limits(aos.eos, args...; kwargs...)
limits(eos::EoSTable, which=Colon()) = begin
    eaxis = ndims(eos.lnEi)==1

    var_min = eaxis ? minimum(eos.lnEi) : minimum(eos.lnT)
    var_max = eaxis ? maximum(eos.lnEi) : maximum(eos.lnT)
    
    return if which == Colon()
        var_min, var_max, minimum(eos.lnRho), maximum(eos.lnRho)
    elseif which == 1
        var_min, var_max
    elseif which == 2
        minimum(eos.lnRho), maximum(eos.lnRho)
    else
        var_min, var_max
    end
end

"""
Return Energy axis of the given EoS.
"""
function EnergyAxis(eos::E; axis=nothing) where {E<:RegularEoSTable}
    eaxis = if isnothing(axis)
        ndims(eos.lnEi) == 1 ? :lnEi : :lnT
    else
        axis
    end

    axis_val  = getfield(eos, eaxis)
    axis_dims = ndims(axis_val)
    axis_len  = first(size(axis_val))
    uni_range = range(minimum(axis_val), maximum(axis_val), length=axis_len)

    isuni = is_uniform(axis_val)
    #if !isuni
    #    @warn "The chosen Energy grid is not uniform."
    #end

    #if axis_dims > 1
    #    @warn "The chosen Energy grid is 2D. Lookup function will use ScatteredInterpolation.jl. This may cause significant slowdown."
    #end

    EoSAxis(axis_val, axis_dims, axis_len, eaxis, uni_range, isuni)
end

"""
Return Density axis of the given EoS.
"""
function DensityAxis(eos::E) where {E<:RegularEoSTable}
    axis_val  = getfield(eos, :lnRho)
    axis_dims = ndims(axis_val)
    axis_len  = last(size(axis_val))
    uni_range = range(minimum(axis_val), maximum(axis_val), length=axis_len)

    isuni = is_uniform(axis_val)
    #if !isuni
    #    @warn "The chosen Density grid is not uniform."
    #end

    #if axis_dims > 1
    #    @warn "The chosen Density grid is 2D. Lookup function will use ScatteredInterpolation.jl. This may cause significant slowdown."
    #end

    EoSAxis(axis_val, axis_dims, axis_len, :lnRho, uni_range, isuni)
end

is_uniform(aos::E) where {E<:AxedEoS}         = is_uniform(aos.energy_axes) & is_uniform(aos.density_axes)
is_uniform(eos::E) where {E<:RegularEoSTable} = is_uniform(EnergyAxis(eos)) & is_uniform(DensityAxis(eos))
is_uniform(a::E)   where {E<:EoSAxis}         = a.is_uniform

"""
    is_uniform(array_from_range)

Determines if an array, that has been constructed from a range collect, 
is equally spaced. 
Note: This actually is quite tricky for larger arrays. The ith entry is 
constructed from \n
    Σj a[1] + j⋅step, \n
until i, however this means that the uncertainty
of the step as well as the uncertainty of the first point need to be propageted.
Even worse, when checking for equal steps after creation, when the initial
step (eps precision) is not known anymore, it can only be computed from the 
entries of the array itself, which have eps uncertainty aswell. Ignoring this, 
the error propagates as \n
    Δstep² = eps(step)².\n
The error of the ith entry is then \n
    Δa[i]² = eps(a[1])² + i⋅Δstep² + eps(a[i])².\n
And hence the error of the distance of two array elements minus the step is \n
    tol = Δa[i-1]² + Δa[i]² + Δstep²),\n
which grows when the array gets longer. We require only that half of 
all significant places match, i.e. \n
    |Δa - step| ≤ √tol.\n
This may need to be improved in the future.
"""
is_uniform(values::A) where {N, T, A<:AbstractArray{T, N}} = begin
    if ndims(values) > 1 
        return false
    end

    if all(diff(values) .≈ first(diff(values)))
        return true
    end

    d::T  = T(values[2] .- values[1])
    d2::T = T(0.0)

    ## eps(d) is the numerical precision of the step
    ## when creating the arrays, the step error is 
    ## added up to the initial error of the floating point value
    ## So each value in the array has a higher uncertainty.
    uni  = true
    lv   = length(values)
    told = eps(d)^2 #eps(values[1])^2 + eps(values[2])^2 + eps(d)^2  # error of the difference 

    tol(i) = eps(values[1])^2 + i*told + eps(values[i])^2

    @inbounds for i in 3:lv
        ##           ϵ(i-1)   + ϵ(i)   + ϵ(d)
        tol_i = sqrt(tol(i-1) + tol(i) +told)
        d2    = values[i] - values[i-1]

        if !isapprox(d2, d, atol=sqrt(tol_i))
            uni = false
            break
        end
    end

    uni
end

@inline EnergyAxis(aos::AxedEoS)  = aos.energy_axes
@inline DensityAxis(aos::AxedEoS) = aos.density_axes

"""
    is_energy(aos)

Check if the energy axis of this EoS is internal energy or not.
"""
is_internal_energy(aos::AxedEoS) = EnergyAxis(aos).name == :lnEi
dependent_energy_variable(aos::AxedEoS) = is_internal_energy(aos) ? :lnT : :lnEi
energy_variable(aos::AxedEoS) = EnergyAxis(aos).name

## Convenience
meshgrid(aos::A) where {A<:AxedEoS}     = meshgrid(EnergyAxis(aos).values, DensityAxis(aos).values)
meshgrid(aos::A) where {A<:RegularEoSTable} = meshgrid(@axed(aos))
size(aos::A) where {A<:AxedEoS}         = (EnergyAxis(aos).length, DensityAxis(aos).length)
size(eos::A) where {A<:RegularEoSTable} = size(@axed(eos))







#= Lookup functions =#

"""
    lookup(eos, what::Symbol, rho, var, args...)
Lookup "what" in the EoS, return the value spliced as args... indicate.
"""
lookup(eos::E, what::Symbol, rho, var, args...) where {E1<:EoSTable, E2<:AxedEoS, E<:Union{E1, E2}} = lookup(lookup_function(eos, what, args...), rho, var)
lookup(eos::E, opacities::OpacityTable, what::Symbol, rho, var, args...)  where {E1<:EoSTable, E2<:AxedEoS, E<:Union{E1, E2}} = begin
    wavelength_window = if (length(args) > 0)
        args[1]
    else
        1:length(opacities.λ)
    end

    if (length(wavelength_window)==1) | (ndims(getfield(opacities, what))==2)
        lookup(lookup_function(eos, opacities, what, args...), rho, var)
    else # lookup all of the wavelength points in the window
        v = zeros(eltype(opacities.κ), size(rho)..., length(opacities.λ[wavelength_window]))
        c = 1
        for i in eachindex(opacities.λ)[wavelength_window]
            setindex!(v, lookup(lookup_function(eos, opacities, what, i), rho, var), (Base.Colon() for j in 1:ndims(rho))..., c)
            c += 1
        end

        v
    end
end

## Generic lookup procedure
@inline lookup(lf::LookupFunction, rho::V, var::V) where {V<:AbstractFloat} = lf.f(var, rho)
@inline lookup(lf::LookupFunction, rho::V, var::V) where {V2<:AbstractFloat, N, V<:AbstractArray{V2,N}} = lf.f.(var, rho)


## Specific methods for specific lookups
@inline lookup(lf::F, rho::V, var::V) where {F<:ScatteredLookup, V<:AbstractFloat} = convert(V, first(lf.f([var], [rho])))
@inline lookup(lf::F, rho::V, var::V) where {F<:ScatteredLookup, V2<:AbstractFloat, N, V<:AbstractArray{V2, N}} = convert.(V2, lf.f(var, rho))


## Axed Lookup function (speed improvements)
lookup_function(eos::EA, what::Symbol, args...) where {EA<:AxedEoS} = begin
    ea = EnergyAxis(eos)
    da = DensityAxis(eos)   

    return if what == ea.name 
        ## Check if the given what is one of the axes
        ## if yes, the lookup function that is created from this should use bisect to lookup!
        other_name = is_internal_energy(eos) ? :lnT : :lnEi
        BisectionLookup((y, x) -> bisect(eos; lnRho=x, (eos.energy_axes.name=>limits(eos.eos, 1), other_name=>y)...))

    elseif is_uniform(eos)
        second_axis = ea.unirange
        rho_axis    = da.unirange
        UniformLookup(linear_interpolation((second_axis, rho_axis), view(getfield(eos.eos, what), :, :, args...), extrapolation_bc=Line()))

    elseif ea.dimension > 1 | da.dimension > 1
        if !scipy_loaded[]
            load_scipy_interpolate!()
        end

        second_axis = if ea.dimension >1 
            ea.values
        else
            puff_up!(zeros(eltype(eos.eos.lnPg), size(eos.eos.lnPg)...), ea.values, 1)
        end

        rho_axis = if da.dimension >1 
            da.values
        else
            puff_up!(zeros(eltype(eos.eos.lnPg), size(eos.eos.lnPg)...), da.values, 2)
        end

        val = view(getfield(eos.eos, what), :, :, args...)

        ## For now this only available using python routines, so we use scipy for this
        ip = (args...) -> pyconvert(
            Array{eltype(eos.eos.lnPg)}, 
            scipy_interpolate.LinearNDInterpolator(stretch(second_axis, rho_axis, val)...)(args...)
        )

        ScatteredLookup(ip)
    else
        second_axis = Interpolations.deduplicate_knots!(copy(ea.values), move_knots=true)
        rho_axis = Interpolations.deduplicate_knots!(copy(da.values), move_knots=true)
        @warn "The EoS seems to be non-uniform! Gridded lookup function will be used."
        GriddedLookup(extrapolate(interpolate((second_axis, rho_axis), view(getfield(eos.eos, what), :, :, args...), Gridded(Linear())), Line()))
    end
end
lookup_function(eos::EA, opacities::O, what::Symbol, args...) where {EA<:AxedEoS, O<:OpacityTable} = begin
    ip = if is_uniform(eos)
        second_axis = EnergyAxis(eos).unirange
        rho_axis    = DensityAxis(eos).unirange
        linear_interpolation((second_axis, rho_axis), log.(view(getfield(opacities, what), :, :, args...)), extrapolation_bc=Line())
    elseif ea.dimension > 1 | da.dimension > 1
        if !scipy_loaded[]
            load_scipy_interpolate!()
        end

        second_axis = if ea.dimension >1 
            ea.values
        else
            puff_up!(zeros(eltype(eos.eos.lnPg), size(eos.eos.lnPg)...), ea.values, 1)
        end

        rho_axis = if da.dimension >1 
            da.values
        else
            puff_up!(zeros(eltype(eos.eos.lnPg), size(eos.eos.lnPg)...), da.values, 2)
        end

        val = log.(view(getfield(opacities, what), :, :, args...))

        ## For now this only available using python routines, so we use scipy for this
        (args...) -> pyconvert(Array{eltype(eos.eos.lnPg)},
            scipy_interpolate.LinearNDInterpolator(stretch(second_axis, rho_axis, val)...)(args...)
        )

    else
        second_axis = Interpolations.deduplicate_knots!(copy(EnergyAxis(eos).values), move_knots=true)
        rho_axis    = Interpolations.deduplicate_knots!(copy(DensityAxis(eos).values), move_knots=true)
        extrapolate(interpolate((second_axis, rho_axis), log.(view(getfield(opacities, what), :, :, args...)), Gridded(Linear())), Line())
    end

    if is_uniform(eos) 
        UniformLookup((x1, x2) -> exp(ip(x1, x2)))
    elseif ea.dimension > 1 | da.dimension > 1
        ScatteredLookup((x1, x2) -> exp.(ip(x1, x2)))
    else
        GriddedLookup((x1, x2) -> exp(ip(x1, x2)))
    end
end


## Regular Lookup function
lookup_function(eos::EB, what::Symbol, args...) where {EB<:RegularEoSTable} = begin
    lookup_function(@axed(eos), what, args...)
    #GriddedLookup(extrapolate(interpolate((EnergyAxis(eos).values, eos.lnRho), view(getfield(eos, what), :, :, args...), Gridded(Linear())), Line()))
end
lookup_function(eos::EB, opacities::O, what::Symbol, args...) where {EB<:RegularEoSTable, O<:OpacityTable} = begin
    lookup_function(@axed(eos), opacities, what, args...)
    #ip = extrapolate(interpolate((EnergyAxis(eos).values, eos.lnRho), log.(view(getfield(opacities, what), :, :, args...)), Gridded(Linear())), Line())
    #GriddedLookup((x1, x2) -> exp(ip(x1, x2)))
end


## Specific lookup function (to bypass form requirement)
lookup_function(::Type{ScatteredLookup}, aos::A, what::Symbol) where {A<:AxedEoS} = begin
    if !scipy_loaded[]
        load_scipy_interpolate!()
    end

    ScatteredLookup(
        (args...)->pyconvert(Array{typeof(aos.eos.lnPg)},
            scipy_interpolate.LinearNDInterpolator(stretch(EnergyAxis(aos).values, DensityAxis(aos).values), getfield(aos.eos, which))(args...)
            )
        )
end 
lookup_function(::Type{ScatteredLookup}, aos::A, opa::B, what::Symbol, args...) where {A<:AxedEoS, B<:OpacityTable} = begin
    if !scipy_loaded[]
        load_scipy_interpolate!()
    end

    ip = ScatteredLookup(scipy_interpolate.LinearNDInterpolator(stretch(EnergyAxis(aos).values, DensityAxis(aos).values), log.(view(getfield(opa, which), :, i))))
    (x,y)->exp.(pyconvert(Array{typeof(aos.eos.lnPg)}, ip(x,y)))
end


## Save interface (maybe a bit slower though)
lookup(eos::A, what::Symbol; lnRho, kwargs...) where {A<:RegularEoSTable} = lookup(@axed(eos), what; lnRho=lnRho, kwargs...)
lookup(eos::A, opa::B, what::Symbol; iλ, lnRho, kwargs...) where {A<:RegularEoSTable, B<:OpacityTable} = lookup(@axed(eos), opa, what; iλ=iλ, lnRho=lnRho, kwargs...)

function lookup(aos::A, what::Symbol; lnRho, kwargs...) where {A<:AxedEoS}
    ename   = energy_variable(aos)
    givenE  = first(keys(kwargs))
    givenEV = kwargs[givenE]
    otherE  = givenE == :lnT ? :lnEi : :lnT

    if (givenE != ename) & (what != ename)
        error("The given E variable is not the Eaxis of the table. If you want to bisect, please pass the E variable as second argument.")
    end

    if (givenE != dependent_energy_variable(aos)) & (what == ename)
        error("The variable to search is the E axis of the table. Please pass the complementary E variable to bisect.")
    end

    lookup(aos, what, lnRho, givenEV)
end
function lookup(aos::A, opa::B, what::Symbol; iλ, lnRho, kwargs...) where {A<:AxedEoS, B<:OpacityTable}
    ename   = energy_variable(aos)
    givenE  = first(keys(kwargs))
    givenEV = kwargs[givenE]
    otherE  = givenE == :lnT ? :lnEi : :lnT

    if (givenE != ename) & (what != ename) & (givenE != dependent_energy_variable(aos))
        error("The given E variable is not the Eaxis of the table. If you want to bisect, please pass the E variable as second argument.")
    end

    givenEV = if (givenE != ename) & (what != ename)
        lookup(aos, ename, lnRho, givenEV)
    else
        givenEV
    end

    iλ > length(opa.λ) && error("Lambda out of bounds for opacities")

    lookup(aos, opa, what, lnRho, givenEV, iλ)
end

ZeroLookup() = ZeroLookup((args...; kwargs...)->0.0)








#= Looking up axes values from others =#

"""
Bisect the EOS to get a value of either of the variables 
on the basis of the other one an one parameter. Enter the names and values as kwargs.
The one that shall be found needs to be given a (start, end).
"""
function bisect(aos::E, iterations=50, antilinear=false; kwargs...) where {E<:AxedEoS}
    eaxis = is_internal_energy(aos)
    bisect_var = eaxis ? :lnEi : :lnT
    given_var  = eaxis ? :lnT  : :lnEi

    # now loop and find the requested values
    bisect_res = 0
    b0, b1     = kwargs[bisect_var]
    bnds       = [b0, b1]
   
    Tp   = eltype(kwargs[:lnRho]) 
    args = Tp[0.0,0.0]
    args[1] = kwargs[given_var]

    f = lookup_function(aos, given_var)


    j = 0
    l = lookup(f, Tp(kwargs[:lnRho]), Tp(first(bnds)))
    while l > kwargs[given_var]
        bnds[1] -= 0.01*(bnds[2]-bnds[1])

        l = lookup(f, Tp(kwargs[:lnRho]), Tp(first(bnds)))
        j+=1
        j==20 && break
    end

    j = 0
    while lookup(f, Tp(kwargs[:lnRho]), Tp(last(bnds))) < kwargs[given_var]
        bnds[2] += 0.01*(bnds[2]-bnds[1])
        j+=1

        j==20 && break
    end

    b0, b1 = Tp.(bnds)



    for i in 1:iterations
        bisect_res = (b0+b1)/2
        args[2] = bisect_res

        para_res = lookup(f, kwargs[:lnRho], args[2])
        para_res = size(para_res) == 0 ? para_res : para_res[1]
     
        if !(antilinear)
            if para_res > args[1]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        else
            if para_res < args[1]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        end
    end

    return bisect_res
end

function bisect(lnρ::AbstractVector, lnT::AbstractVector, aos::E) where {E<:AxedEoS}
    lnEi = similar(lnT)
    emin,emax = limits(aos.eos, 1)

    var = aos.energy_axes.name
    othervar = var == :lnEi ? :lnT : :lnEi

    for i in eachindex(lnρ)
        lnEi[i] = bisect(aos, lnRho=lnρ[i], (othervar=>lnT[i], var=>[emin, emax])...)
    end

    #hcat(z, exp.(lnEi), exp.(lnρ));
    lnEi
end

bisect(eos::E, args...; kwargs...) where {E<:RegularEoSTable} = bisect(AxedEoS(eos), args...; kwargs...)
bisect(aos::E, lnRho::F, between=limits(aos.eos, 1); kwargs...) where {E<:AxedEoS, F<:AbstractFloat} = bisect(aos; lnRho=lnRho, (aos.energy_axes.name=>between,)..., kwargs...)






#= Image Filtering for smoothing the EoS =#

"""
    smooth(eos::AxedEoS, radius=30)

Apply a gaussian 1D filter to the energy dimension of the given EoS table.
"""
function smooth(aos::AxedEoS, radius=30)
	kernel = ImageFiltering.Kernel.gaussian((radius,))

    eos = aos.eos
    V = getfield(aos.eos, dependent_energy_variable(aos))

    E_conv = similar(V)
	P_conv = similar(V)
	N_conv = similar(V)
	R_conv = similar(V)
	
	for i in axes(E_conv, 2)
		x=V[:, i]
		E_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnPg[:, i]
		P_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnNe[:, i]
		N_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnRoss[:, i]
		R_conv[:, i] .= ImageFiltering.imfilter(x, kernel)
	end

    if is_internal_energy(aos)
	    SqEoS(eos.lnRho, E_conv, eos.lnEi, P_conv, R_conv, N_conv)
    else
        SqEoS(eos.lnRho, eos.lnT, E_conv, P_conv, R_conv, N_conv)
    end
end

"""
    smooth(opa::SqOpacity, radius=30)

Apply a gaussian 1D filter to the energy dimension of the given opacity table.
"""
function smooth(opa::SqOpacity, radius=30)
	kernel = ImageFiltering.Kernel.gaussian((radius,))

    k_conv = similar(opa.κ)
	kr_conv = similar(opa.κ_ross)
	src_conv = similar(opa.κ)
	
	for l in axes(k_conv, 3)
		for i in axes(k_conv, 2)
			x=log.(opa.κ[:, i, l])
			k_conv[:, i, l] .= ImageFiltering.imfilter(x, kernel)
			
			x=log.(opa.src[:, i, l])
			src_conv[:, i, l] .= ImageFiltering.imfilter(x, kernel)
		end
	end
	for i in axes(k_conv, 2)
		x=log.(opa.κ_ross[:, i])
		kr_conv[:, i] .= ImageFiltering.imfilter(x, kernel)
	end
	
	SqOpacity(exp.(k_conv), exp.(kr_conv), exp.(src_conv), opa.λ, opa.optical_depth)
end

"""
    smooth(eos::SqEoS, opa::SqOpacity; eos_radius=30, opa_radius=30)

Apply a gaussian 1D filter to the energy dimension of the given EoS and opacity table.
"""
smooth(eos::SqEoS, opa::SqOpacity; eos_radius=30, opa_radius=30) = smooth(@axed(eos), eos_radius), smooth(opa, opa_radius)
smooth(eos::AxedEoS, opa::SqOpacity; eos_radius=30, opa_radius=30) = smooth(eos, eos_radius), smooth(opa, opa_radius)






#= Writing in Dispatch format =#

function tabparam(folder::String, EOSTableFile, RhoEiRadTableFile, tabparam_file, nEiBin, nRhoBin, nRadBins, EiMin, EiMax, RhoMin, RhoMax)
    content="&TABPAR\n"*
            "EOSTableFile = '$(EOSTableFile)'\n"*
            "RhoEiRadTableFile = '$(RhoEiRadTableFile)'\n"*
            "NeTgRadTableFile = 'netg_radtab.dat'\n"*
            "nRhoBin = $(nRhoBin)\n"*
            @sprintf("RhoMin = %.7E\n", RhoMin)*
            @sprintf("RhoMax = %.7E\n", RhoMax)*
            "nEiBin = $(nEiBin)\n"*
            @sprintf("EiMin =  %.7E\n", EiMin)*
            @sprintf("EiMax =  %.7E\n", EiMax)*
            "nRadBins = $(nRadBins)\n"*
            "/\n"

    open(joinpath(folder, tabparam_file), "w") do f
        write(f, content)
    end

    nothing
end

tabparam(eos::EoSTable, nradbins, folder::String; eos_file="eostable.dat", opacity_file="rhoei_radtab.dat", tabparam_file="tabparam.in") = tabparam(folder, eos_file, opacity_file, tabparam_file, size(eos.lnPg)..., nradbins, exp.(limits(eos))...)

"""
Write the EoS + binned opacities in the same format as in Tabgen, so that it can be read by dispatch.
"""
function for_dispatch(eos::EoSTable, χ, S, ϵ, folder::String; name="")
    !isdir(folder) && mkdir(folder) 

    f = FortranFile(joinpath(folder, "rhoei_radtab$(name).dat"), "w", access="direct", recl=prod(size(S))*4)
    FortranFiles.write(f, rec=1, log.(ϵ))
    FortranFiles.write(f, rec=2, log.(S))
    FortranFiles.write(f, rec=3, log.(χ))
    close(f)

    eos_table = zeros(eltype(eos.lnT), size(eos.lnT)..., 4)
    eos_table[:, :, 1] = eos.lnPg
    eos_table[:, :, 2] = exp.(eos.lnT)
    eos_table[:, :, 3] = eos.lnNe
    eos_table[:, :, 4] = eos.lnRoss

    f = FortranFile(joinpath(folder, "eostable$(name).dat"), "w", access="direct", recl=prod(size(eos.lnT))*4*4)
    FortranFiles.write(f, rec=1, eos_table)
    close(f)

    tabparam(eos, size(χ, 3), folder, eos_file="eostable$(name).dat", opacity_file="rhoei_radtab$(name).dat", tabparam_file="tabparam$(name).in")
end

for_dispatch(eos::EoSTable, opacities::OpacityTable, folder::String; name="") = begin
    eos_table_name = folder

    for_dispatch(eos, opacities.κ, opacities.src, opacities.κ_ross, folder)
    
    save(opacities, joinpath(eos_table_name, "binned_opacities$(name).hdf5"))
    save(eos, joinpath(eos_table_name, "eos$(name).hdf5"))
end

for_dispatch(aos::AxedEoS, args...; kwargs...) = begin
    if !is_internal_energy(aos)
        error("Energy axis has to be internal energy for DISPATCH. The current axis is $(aos.energy_axes.name)")
    end

    if !aos.energy_axes.is_uniform 
        @warn "A non-uniform Energy axes will be written to DISPATCH. This may cause issues later."
    end

    if !aos.density_axes.is_uniform 
        @warn "A non-uniform Density axes will be written to DISPATCH. This may cause issues later."
    end

    for_dispatch(aos.eos, args...; kwargs...)
end

function save_tables(aos, opa, eos_table_name, table_folder=nothing)
    # Save everything in the dispatch format
    for_dispatch(aos, opa)
    save(opa.opacities, "binned_opacities.hdf5")

    # Move files to the final folder for dispatch
    !isdir(eos_table_name) && mkdir(eos_table_name) 
    mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
    mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
    mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)

    # Copy the eos for convenience. Usually not a big deal because rather small
    if !isnothing(table_folder)
        cp(joinpath(table_folder, "combined_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true);
    end;
end




#= Manipulating opacity tables =#

scale!(opa::SqOpacity; fields...) = begin
    for (fieldname, factor) in fields
        getfield(opa, fieldname) .*= factor
    end
end