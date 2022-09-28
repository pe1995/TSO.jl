abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end

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

OpacityTable = Union{RegularOpacityTable, IrregularOpacityTable} 
EoSTable     = Union{RegularOpacityTable, IrregularEoSTable} 

### Convenience Constructor functions
Opacity(args...; regular=true, kwargs...) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...; regular=true, kwargs...)     = regular ? RegularEoSTable(args...; kwargs...)     : IrregularEoSTable(args...; kwargs...)

RegularOpacityTable(κ::Array, κ_ross::Array, src::Array; optical_depth=false) = RegularOpacityTable(κ, κ_ross, src, optical_depth) 