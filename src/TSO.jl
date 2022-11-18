module TSO

# Imports
using BSplineKit
using Interpolations
using Printf
using DelimitedFiles
using FortranFiles
using PyCall
using Glob
using HDF5
using Statistics

# Interface
abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end
abstract type OpacityBins end

# Source files
include("_aux.jl")
include("_communication.jl")
include("_tables.jl")
include("_eos.jl")
include("_opacities.jl")
include("_smooth.jl")
include("_binning.jl")
include("_transfer.jl")

end
