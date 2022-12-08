"""Table containing binned Opacities + Source functions on a regular grid."""
struct RegularOpacityTable{T<:AbstractFloat, 
                            NOp, NRoss, NSrc} <: AbstractRegularTable
    κ             ::Array{T, NOp}
    κ_ross        ::Array{T, NRoss}
    src           ::Array{T, NSrc}
    λ             ::Array{T, 1}
    optical_depth ::Bool
end

"""Table containing EoS properties on a regular grid."""
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

"""Table containing binned Opacities + Source functions on an irregular grid."""
struct IrregularOpacityTable <: AbstractIrregularTable
end

"""Table containing EoS properties on an irregular grid."""
struct IrregularEoSTable <: AbstractIrregularTable
end


### Aliases
OpacityTable = Union{RegularOpacityTable, IrregularOpacityTable} 
EoSTable     = Union{RegularEoSTable, IrregularEoSTable} 

SqOpacity = RegularOpacityTable
SqEoS     = RegularEoSTable



### Convenience Constructor functions
Opacity(args...; regular=true, kwargs...) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...; regular=true, kwargs...)     = regular ? RegularEoSTable(args...; kwargs...)     : IrregularEoSTable(args...; kwargs...)

RegularOpacityTable(κ::Array, κ_ross::Array, src::Array; optical_depth=false) = RegularOpacityTable(κ, κ_ross, src, optical_depth) 



### General functions
limits(eos::EoSTable) = begin
    eaxis = ndims(eos.lnT)>1

    var_min = eaxis ? minimum(eos.lnEi) : minimum(eos.lnT)
    var_max = eaxis ? maximum(eos.lnEi) : maximum(eos.lnT)
    
    var_min, var_max, minimum(eos.lnRho), maximum(eos.lnRho)
end



### Writing in Dispatch format
function tabparam(EOSTableFile, RhoEiRadTableFile, nEiBin, nRhoBin, nRadBins, EiMin, EiMax, RhoMin, RhoMax)
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
            "/"

    open("tabparam.in", "w") do f
        write(f, content)
    end

    nothing
end

tabparam(eos::EoSTable, nradbins; eos_file="eostable.dat", opacity_file="rhoei_radtab.dat") = tabparam(eos_file, opacity_file, size(eos.lnPg)..., nradbins, exp.(TSO.limits(eos))...)

"""
Write the EoS + binned opacities in the same format as in Tabgen, so that it can be read by dispatch.
"""
function for_dispatch(eos::EoSTable, χ, S, ϵ)
    f = FortranFile("rhoei_radtab.dat", "w", access="direct", recl=prod(size(S))*4)
    FortranFiles.write(f, rec=1, ϵ)
    FortranFiles.write(f, rec=2, S)
    FortranFiles.write(f, rec=3, χ)
    close(f)

    eos_table = zeros(eltype(eos.lnT), size(eos.lnT)..., 4)
    eos_table[:, :, 1] = eos.lnPg
    eos_table[:, :, 2] = exp.(eos.lnT)
    eos_table[:, :, 3] = eos.lnNe
    eos_table[:, :, 4] = eos.lnRoss

    f = FortranFile("eostable.dat", "w", access="direct", recl=prod(size(eos.lnT))*4*4)
    FortranFiles.write(f, rec=1, eos_table)
    close(f)

    tabparam(eos, size(χ, 3))
end

for_dispatch(eos::EoSTable, opacities::OpacityTable) = for_dispatch(eos, opacities.κ, opacities.src, opacities.κ_ross)
for_dispatch(eos::EoSTable, opacities::OpacityTable, folder::String) = begin
    for_dispatch(eos, opacities)

    save(opacities, "binned_opacities.hdf5")
    save(eos, "binned_eos.hdf5")

    # Move files to the final folder for dispatch
    eos_table_name = folder
    !isdir(eos_table_name) && mkdir(eos_table_name) 

    mv("tabparam.in",           joinpath(eos_table_name, "tabparam.in"),           force=true)
    mv("eostable.dat",          joinpath(eos_table_name, "eostable.dat"),          force=true)
    mv("rhoei_radtab.dat",      joinpath(eos_table_name, "rhoei_radtab.dat"),      force=true)
    mv("binned_opacities.hdf5", joinpath(eos_table_name, "binned_opacities.hdf5"), force=true)
    mv("binned_eos.hdf5",       joinpath(eos_table_name, "eos.hdf5"),              force=true);
end