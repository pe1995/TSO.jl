#================================================
======= AESOPUS 2.0 Opacity + EoS ===============
================================================#

## Constructors

"""
    AesopusEoS(path)

Read an Aesopus EoS from a HDF5 file. 
Reads it directly into the squaregas format.
"""
function AesopusEoS(path)
    f = HDF5.h5open(path, "r")

    lnT    = log.(HDF5.read(f["T"]))
    lnEi   = HDF5.read(f["lnEi"])
    lnRoss = HDF5.read(f["lnKross"])
    lnNe   = HDF5.read(f["lnNe"])
    lnPg   = HDF5.read(f["lnPg"])
    lnRho  = HDF5.read(f["lnRho"])

    close(f)

    @assert all(size(lnPg) .== (length(lnT), length(lnRho)))

    RegularEoSTable(lnRho, lnT, lnEi, lnPg, lnRoss, lnNe)
end