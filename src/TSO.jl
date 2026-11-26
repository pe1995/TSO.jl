module TSO

# ============================================================================
# Imports
# ============================================================================
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

# ============================================================================
# Abstract Interface
# ============================================================================
abstract type AbstractModel end

# ============================================================================
# User interface
# ============================================================================
export EquationOfState, availableGrids, migrateTable, scattering, units
export sample, sample!, weights, weights!
export TableInfo, find_table_info, show_composition
export binning, for_dispatch, for_dispatch_T
export rosseland_optical_depth, optical_depth

# ============================================================================
# Python libs
# ============================================================================
const scipy_interpolate = PythonCall.pynew()
const scipy_loaded = Ref(false)
const numpy = PythonCall.pynew()
const numpy_loaded = Ref(false)
const generalTimer = TimerOutput()

# ============================================================================
# Source files
# ============================================================================
include("_aux.jl")
include("_eos_info.jl")
include("_equation_of_state.jl")
include("_tables.jl")
include("_sampling.jl")
include("_planck.jl")
include("_models.jl")
include("_eos.jl")
include("_binning.jl")
include("_adiabat.jl")
include("_from_m3d.jl")
include("_temp2eint.jl")
include("_extend.jl")
include("_dispatch_hooks.jl")

# ============================================================================
# List of timers
# ============================================================================
const timers = [
    rosseland_time,
    optical_depth_time,
    formation_height_time
]

# ============================================================================
# Deprecations
# ============================================================================

end
