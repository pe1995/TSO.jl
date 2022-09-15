abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end

"""Table containing binned Opacities + Source functions on a regular grid."""
struct RegularOpacityTable <: AbstractRegularTable
end

"""Table containing EoS properties on a regular grid."""
struct RegularEoSTable{F<:AbstractFloat} <: AbstractRegularTable
    lnRho  :: Array{F, 3}
    lnT    :: Array{F, 3}
    lnEi   :: Array{F, 3}
    lnRoss :: Array{F, 3}
    lnNe   :: Array{F, 3}
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