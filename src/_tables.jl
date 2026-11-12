struct RegularOpacityTable{T<:AbstractFloat,
                            NOp, NRoss, NSrc} <: AbstractRegularTable
    κ             ::Array{T, NOp}
    κ_ross        ::Array{T, NRoss}
    src           ::Array{T, NSrc}
    λ             ::Array{T, 1}
    optical_depth ::Bool
end

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

struct IrregularOpacityTable <: AbstractIrregularTable
end

struct IrregularEoSTable <: AbstractIrregularTable
end

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

# ============================================================================
# Aliases
# ============================================================================

OpacityTable = Union{RegularOpacityTable, IrregularOpacityTable}
EoSTable     = Union{RegularEoSTable, IrregularEoSTable}
SqOpacity    = RegularOpacityTable
SqEoS        = RegularEoSTable

# ============================================================================
# Convenience Constructor functions
# ============================================================================

Opacity(args...; regular=true, kwargs...) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...;     regular=true, kwargs...) = regular ? RegularEoSTable(args...;     kwargs...) : IrregularEoSTable(args...;     kwargs...)

RegularOpacityTable(κ::Array, κ_ross::Array, src::Array; optical_depth=false) = RegularOpacityTable(κ, κ_ross, src, optical_depth)

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

# ============================================================================
# General functions
# ============================================================================

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

    EoSAxis(axis_val, axis_dims, axis_len, :lnRho, uni_range, isuni)
end

is_uniform(aos::E) where {E<:AxedEoS}         = is_uniform(aos.energy_axes) & is_uniform(aos.density_axes)
is_uniform(eos::E) where {E<:RegularEoSTable} = is_uniform(EnergyAxis(eos)) & is_uniform(DensityAxis(eos))
is_uniform(a::E)   where {E<:EoSAxis}         = a.is_uniform

is_uniform(values::A) where {N, T, A<:AbstractArray{T, N}} = begin
    if ndims(values) > 1
        return false
    end

    if all(diff(values) .≈ first(diff(values)))
        return true
    end

    d::T  = T(values[2] .- values[1])
    d2::T = T(0.0)

    uni  = true
    lv   = length(values)
    told = eps(d)^2

    tol(i) = eps(values[1])^2 + i*told + eps(values[i])^2

    @inbounds for i in 3:lv
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

# ============================================================================
# Convenience
# ============================================================================
meshgrid(aos::A) where {A<:AxedEoS}     = meshgrid(EnergyAxis(aos).values, DensityAxis(aos).values)
meshgrid(aos::A) where {A<:RegularEoSTable} = meshgrid(@axed(aos))
size(aos::A) where {A<:AxedEoS}         = (EnergyAxis(aos).length, DensityAxis(aos).length)
size(eos::A) where {A<:RegularEoSTable} = size(@axed(eos))

# ============================================================================
# Image Filtering for smoothing the EoS
# ============================================================================

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

smooth(eos::SqEoS, opa::SqOpacity; eos_radius=30, opa_radius=30) = smooth(@axed(eos), eos_radius), smooth(opa, opa_radius)
smooth(eos::AxedEoS, opa::SqOpacity; eos_radius=30, opa_radius=30) = smooth(eos, eos_radius), smooth(opa, opa_radius)

# ============================================================================
# Writing in Dispatch format
# ============================================================================

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

function dispatch_axes(folder::String, lnRho, lnE; name="")
    open(joinpath(folder, "eosaxes$(name).dat"), "w") do f
        write(f, Int32(length(lnRho)), Int32(length(lnE)))
        write(f, Float64.(collect(lnRho)))
        write(f, Float64.(collect(lnE)))
    end
end

"""
Write the EoS + binned opacities in the same format as in Tabgen, so that it can be read by dispatch.
"""
function for_dispatch(eos::EoSTable, χ::AbstractArray, S::AbstractArray, ϵ::AbstractArray, folder::String; name="")
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

    dispatch_axes(folder, eos.lnRho, eos.lnEi, name=name)

    tabparam(eos, size(χ, 3), folder, eos_file="eostable$(name).dat", opacity_file="rhoei_radtab$(name).dat", tabparam_file="tabparam$(name).in")
end

"""
Write the NATIVE rho-T binned opacity (κ) and source (S) in DISPATCH format, for the
`square_gas_T` EOS variant. Companion to `for_dispatch`, but on the temperature grid:
only the opacity/source radiation table (`rhoei_radtab_T.dat`) and its `tabparam_T.in`
are written — NO `eostable_T.dat`, because `square_gas_T` obtains all thermodynamics
(T, pg, ne, ross, gamma) from the standard rho-E `eostable.dat` and looks up κ/S on this
rho-T table using that T. `eos` here is the native rho-T table (`aos.eos`): `eos.lnT`
is the 1-D temperature axis and `size(eos.lnPg) == (nT, nRho)`, so `tabparam` auto-fills
the "Ei" namelist slots with the log-T extrema (nEiBin=nT, EiMin/EiMax=Tmin/Tmax).
"""
function for_dispatch_T(eos::EoSTable, χ, S, ϵ, folder::String; name="_T")
    !isdir(folder) && mkdir(folder)

    eax = EnergyAxis(eos)
    if !eax.is_uniform
        error("for_dispatch_T: temperature axis ($(eax.name)) is not uniform in log; " *
              "the DISPATCH square_gas_T reader assumes a uniform log-T grid.")
    end

    f = FortranFile(joinpath(folder, "rhoei_radtab$(name).dat"), "w", access="direct", recl=prod(size(S))*4)
    FortranFiles.write(f, rec=1, log.(ϵ))
    FortranFiles.write(f, rec=2, log.(S))
    FortranFiles.write(f, rec=3, log.(χ))
    close(f)

    dispatch_axes(folder, eos.lnRho, eos.lnT, name=name)

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
    for_dispatch(aos, opa)
    save(opa.opacities, "binned_opacities.hdf5")

    !isdir(eos_table_name) && mkdir(eos_table_name)
    mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
    mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
    isfile("eosaxes.dat") && mv("eosaxes.dat", joinpath(eos_table_name, "eosaxes.dat"), force=true)
    mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)

    if !isnothing(table_folder)
        cp(joinpath(table_folder, "combined_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true);
    end;
end

# ============================================================================
# Manipulating opacity tables
# ============================================================================

scale!(opa::SqOpacity; fields...) = begin
    for (fieldname, factor) in fields
        getfield(opa, fieldname) .*= factor
    end
end
