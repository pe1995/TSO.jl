module TSO

#= Imports =#
import BSplineKit
using Interpolations
using Printf
using DelimitedFiles
using FortranFiles
using PythonCall
using Glob
using HDF5
using Statistics
using Mmap
using Clustering
using ThreadsX
using TimerOutputs
using ImageFiltering
using LoopVectorization
import Clustering.assignments

import Base.size
import Base.broadcastable
import Base.length


#= Abstract Interface =#
abstract type AbstractTable end
abstract type AbstractRegularTable <:AbstractTable end
abstract type AbstractIrregularTable <:AbstractTable end
abstract type AbstractBinnedOpacities end
abstract type OpacityBins end
abstract type AbstractModel end
abstract type LookupFunction end
abstract type RadiativeTransferSolver end
abstract type AbstractRadiationField end



#= User Interface (TBD) =#
## Tables
export OpacityTable, EoSTable, SqOpacity, SqEoS
export for_dispatch
export @axed, limits, EnergyAxis, DensityAxis, AxedEoS, is_internal_energy

## fast lookup and extended tables interfac
export opacity, wavelength, extended, sample, sample!, weights

## Aux
export save, reload, meshgrid

## Binning 
export TabgenBins, TabgenBinning, StaggerBinning, Co5boldBinning, MURaMBinning, DensityBinning, ClusterBinning
export binning
export ω_midpoint
export lookup, lookup_function, bisect
export tabulate
export reset_bins
export @binned, @binned!

## EoS
export rosseland_opacity!, rosseland_opacity, transfer_rosseland!, upsample
export formation_height, optical_depth

## MARCS
export MARCSOS, MARCSOpacity

## Intpolations
export gridded, uniform, complement, switch_energy, set_limits!

## Turbospectrum
export load_wrapper, load_TS, move_output, import_wrapper, inTS, inWrapper
export babsma!, bsyn!

## Aesopus
export AesopusEoS, AesopusOpacity

## Models
export Model1D, Average3D, @optical

## Transfer
export Solver, Jν, Qr


#= Python libs =#
const scipy_interpolate = PythonCall.pynew()
const scipy_loaded = Ref(false)
const numpy = PythonCall.pynew()
const numpy_loaded = Ref(false)
const generalTimer = TimerOutput()

#= Source files =#
include("_aux.jl")
include("_communication.jl")
include("_tables.jl")
include("_mini_eos.jl")
include("_models.jl")
include("_eos.jl")
include("_smooth.jl")
include("_transfer.jl")
include("_heating.jl")
include("_binning.jl")
include("_legacy_tables.jl")
include("_aesopus.jl")
include("_marcs.jl")
include("_adiabat.jl")
include("_from_m3d.jl")
include("_binning_execution.jl")
include("_extend.jl")
include("_extended_eos.jl")


#= List of timers =#
const timers = [
    rosseland_time, 
    optical_depth_time, 
    formation_height_time,
    binning_time
]

#= Deprecations =#
#@deprecate replaceEoS complement do_binning! box_integrated_v3 box_integrated_v2 box_integrated

end
