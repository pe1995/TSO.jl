module TSO

#= Imports =#
using BSplineKit
using Interpolations
using Printf
using DelimitedFiles
using FortranFiles
using PyCall
using Glob
using HDF5
using Statistics


#= Abstract Interface =#
abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end
abstract type OpacityBins end


#= User Interface (TBD) =#
## Tables
export OpacityTable, EoSTable, SqOpacity, SqEoS
export for_dispatch
export EnergyAxis, DensityAxis

## Turbospectrum
export load_wrapper, load_TS, move_output, import_wrapper, inTS, inWrapper
export babsma!, bsyn!

## Aesopus
export AesopusEoS, AesopusOpacity

## Aux
export save, reload

## Binning 
export TabgenBins, TabgenBinning, StaggerBinning, Co5boldBinning, MURaMBinning
export binning
export ω_midpoint
export lookup, bisect
export tabulate

## EoS
export rosseland_opacity!, rosseland_opacity, unify

## MARCS
export MARCSOS, square, complement


#= Source files =#
include("_aux.jl")
include("_communication.jl")
include("_tables.jl")
include("_eos.jl")
include("_opacities.jl")
include("_smooth.jl")
include("_binning.jl")
include("_transfer.jl")
include("_legacy_tables.jl")
include("_aesopus.jl")
include("_marcs.jl")

#= Deprecations =#
@deprecate replaceEoS complement


end
