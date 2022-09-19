abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end

"""Table containing binned Opacities + Source functions on a regular grid."""
struct RegularOpacityTable <: AbstractRegularTable
end

"""Table containing EoS properties on a regular grid."""
struct RegularEoSTable{F<:AbstractFloat, 
                        NRho<:Integer, NT<:Integer, NEi<:Integer, 
                        NPg<:Integer, NRoss<:Integer, NNe<:Integer} <: AbstractRegularTable
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

# Convenience Constructor functions
Opacity(args...; kwargs..., regular=true) = regular ? RegularOpacityTable(args...; kwargs...) : IrregularOpacityTable(args...; kwargs...)
EoS(args...; kwargs..., regular=true)     = regular ? RegularEoSTable(args...; kwargs...)     : IrregularEoSTable(args...; kwargs...)