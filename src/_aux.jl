## General Utility functions
"""
Split array arr in nsplits roughly equal junks.
If mask=true the indices will be returned.
"""
function split_similar(arr, nsplits; mask=false)
    nsplits = length(arr) < nsplits ? length(arr) : nsplits
    splits  = div(length(arr), nsplits)
    nrest   = length(arr) % nsplits
    split_sizes = Int[splits for _ in 1:nsplits]
    
    for i in 1:nrest
        split_sizes[i] = split_sizes[i] + 1
    end

    split_masks = []
    last_idx    = 0
    for i in 1:nsplits
        append!(split_masks, [[ (last_idx+1:last_idx+split_sizes[i])... ]])
        last_idx = last_idx+split_sizes[i]
    end

    if mask
        return split_masks
    else
        return [arr[mask] for mask in split_masks]
    end
end

function meshgrid(ax...)
    grids = []
    space = Iterators.product(ax...)
    for i in eachindex(ax)
        append!(grids, [getindex.(space, i)])
    end

    grids
end

"""
create random number between a and b
"""
randrange(a,b,args...) = begin
    xmin,xmax = min(a,b),max(a,b)
    rand(args...) .* (xmax-xmin) .+ xmin
end

## Importing of python modules

function _get_help_py(mod ,dir=dirname(@__FILE__))
	sys   = pyimport("sys")
	dir in sys."path" ? nothing : append!(sys."path",[dir])
    pyimport(String(mod))
end

macro pythonHelp(mod)
    mod_e = esc(mod)
    mod_s = :($mod)
    path  = dirname(@__FILE__)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path))
end

macro pythonHelp(mod, dir)
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end

## Saving as HDF5
function save(s::T, path) where {T<:AbstractTable}
    fid = HDF5.h5open(path, "w")

    for fname in fieldnames(typeof(s))
        add_to_hdf5!(fid, fname, getfield(s, fname))
    end

    close(fid)

    path
end

function reload(s::Type{S}, path::String; mmap=false) where {S}
    fid   = HDF5.h5open(path, "r")
    fvals = Any[]

    for (fname, ftype) in zip(fieldnames(s), fieldtypes(s))
        append!(fvals, [get_from_hdf5(ftype, fid, fname; mmap=mmap)])
    end

    close(fid)

    s(fvals...)
end

add_to_hdf5!(fid, fname, val)       = fid["$(fname)"] = val
add_to_hdf5!(fid, fname, val::Bool) = fid["$(fname)"] = Int(val)

get_from_hdf5(::Type{<:Any}, fid, fname; mmap=false) = mmap ? HDF5.readmmap(fid["$(fname)"]) : HDF5.read(fid["$(fname)"])
get_from_hdf5(::Type{Bool},  fid, fname; mmap=false) = Bool(HDF5.read(fid["$(fname)"]))
    

## Constants (In agreement with Tabgen)
const EiExtra    = 5.0                         # extra eV per H atom to add
const KBoltzmann = 1.380658E-16                # Boltzman's cst. [erg/K]
const CLight     = 2.99792458E+10              # Speed of light [cm/s]
const HPlanck    = 6.6260755E-27               # Planck's constant [erg s]
const ev_to_erg  = 1.60218e-12                 # conversion
const HIonPot    = 13.595                      # hydrogen ioniz. potential [eV]
const twohc2     = 2.0e0 *HPlanck*CLight^2
const hc_k       = HPlanck*CLight/KBoltzmann
const aa_to_cm   = 1.0e-8

## Plot default setup that kind-of works
#=
import PyCall.rc
rc("mathtext", fontset="cm", default="regular")
rc("font", family="monospace", size=12)
rc("text", usetex=false)
rc("xtick", direction="in", top=true)
rc("ytick", direction="in", right=true)
rc("xtick.minor", visible=true)
rc("ytick.minor", visible=true)
=#