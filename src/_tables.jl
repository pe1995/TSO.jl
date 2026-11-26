# ============================================================================
# Table axes
# ============================================================================

struct EoSAxis{F<:AbstractFloat, Naxis, I<:Integer, R<:AbstractRange}
    values     ::Array{F, Naxis}
    dimension  ::I
    length     ::I
    name       ::Symbol
    unirange   ::R
    is_uniform ::Bool
end

is_uniform(a::E) where {E<:EoSAxis} = a.is_uniform

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

tabparam(eos::EquationOfState, nradbins, folder::String; eos_file="eostable.dat", opacity_file="rhoei_radtab.dat", tabparam_file="tabparam.in") = tabparam(folder, eos_file, opacity_file, tabparam_file, size(eos.lnPg)..., nradbins, exp.(limits(eos))...)

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
function for_dispatch(eos::EquationOfState, χ::AbstractArray, S::AbstractArray, ϵ::AbstractArray, folder::String; name="")
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
rho-T table using that T. `eos` here is the native rho-T table: `eos.lnT`
is the 1-D temperature axis and `size(eos.lnPg) == (nT, nRho)`, so `tabparam` auto-fills
the "Ei" namelist slots with the log-T extrema (nEiBin=nT, EiMin/EiMax=Tmin/Tmax).
"""
function for_dispatch_T(eos::EquationOfState, χ, S, ϵ, folder::String; name="_T")
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

"""
    for_dispatch(tab, folder; name="")

Write a table on the internal energy grid with its binned opacities for DISPATCH.
"""
for_dispatch(tab::EquationOfState, folder::String; name="") = begin
    is_internal_energy(tab) || error("Energy axis has to be internal energy for DISPATCH. The current axis is $(energy_variable(tab))")
    is_uniform(tab.energy_axis) || @warn "A non-uniform Energy axes will be written to DISPATCH. This may cause issues later."
    is_uniform(tab.density_axis) || @warn "A non-uniform Density axes will be written to DISPATCH. This may cause issues later."
    for_dispatch(tab, tab[:κ], tab[:src], ones(eltype(tab[:src]), size(tab[:src])...), folder; name=name)
end
