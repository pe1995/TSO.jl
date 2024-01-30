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

## Aux
export save, reload, meshgrid

## Binning 
export TabgenBins, TabgenBinning, StaggerBinning, Co5boldBinning, MURaMBinning, DensityBinning, ClusterBinning
export binning
export ω_midpoint
export lookup, lookup_function, bisect
export tabulate
export reset_bins
export @binned

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
const scipy_loaded      = Ref(false)


#= Source files =#
include("_aux.jl")
include("_communication.jl")
include("_tables.jl")
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



#= Deprecations =#
@deprecate replaceEoS complement
#@deprecate regrid replace_axis



#= TODO =#
#=
    MOST IMPORTANT PROBLEMS:
    ------------------------
(A) The resolution when switching from T
    to E grids needs to be increased by a lot
    because otherwise horizontal parts in the T-E
    plane are very poorly interpolated!

    -> A possible solution is to keep the tables 
    that are on the T grid on there as long as possible!
    When they are binned the tables are tiny anyways, 
    so then a interpolation to a very dense E grid is possible.

    -> For this one has to make sure that every function in here is
    capable of working with a temperature grid, i.e. does not
    require internal energy as before! Use the AxedEoS everywhere
    to check if the grid is T or E.

(B) All functions work only if the input tables have some sort
    of symmetry, i.e. are gridded in at least one dimension.
    It would be much nicer if the lookup infrastructure would
    be able to interpolate on any grid (also scattered),
    so that every table however it looks like can be interpolated
    to any grid by one sinple function. The EoS just needs to be 
    Axed, then the lookup knows how to interpolate, so that you can
    just loop though whatever new grid and lookup all the quantites for it.
    There should be function that accepts any new axes and a function that
    just makes it regular in the alredy existing axes.

    UPDATE 1: Scattered Interpolation is taken over from scipy for now
              -> there is no extrapolation yet, howver can be added to
              the square function (tbd, not there yet) by extrapolating 
              constant rows
=#


end
