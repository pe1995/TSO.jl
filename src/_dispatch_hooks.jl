# =============================================================================
# EoS hooks for the MUST.jl snapshot conversion
# =============================================================================
#
# `MUST.snapshotBox` takes an
# `eos_reader` + `lookup_generator` pair to use an EoS from outside instead:
#
#   MUST.snapshotBox(
#       snap; folder=rundir,
#       eos_reader=TSO.dispatch_eos_reader(eos_dir),
#       lookup_generator=TSO.dispatch_lookup_generator
#   )
#
# where `eos_dir` is the table directory the run's `table_loc` points at.
#
# =============================================================================

const dispatch_eos_quantities = [:T, :kr, :Pg, :Ne]
const dispatch_eos_variables = Dict(:T=>:lnT, :kr=>:lnRoss, :Pg=>:lnPg, :Ne=>:lnNe)

"""
    dispatch_eos_reader(eos_path; quantities=dispatch_eos_quantities)

`eos_reader` for `MUST.snapshotBox`, loading the energy based TSO EoS at
`eos_path` (a table directory or the `eos.hdf5` itself).
"""
dispatch_eos_reader(eos_path; quantities=dispatch_eos_quantities) = (run) -> begin
    path = isdir(eos_path) ? joinpath(eos_path, "eos.hdf5") : eos_path
    isfile(path) || error("No TSO EoS at $(path).")

    eos = extended(reload(SqEoS, path))
    energy_variable(eos.eos) == :lnEi || error(
        "$(path) is not an energy based table, the conversion looks up in (lnRho, lnEi)."
    )

    eos, quantities
end

"""
    dispatch_lookup_generator(eos, q)

`lookup_generator` for `MUST.snapshotBox`. Returns the payload for quantity `q`
and the shared kernel that fills it in.
"""
function dispatch_lookup_generator(eos::ExtendedEoS, q)
    coefs = Ref{Any}(nothing)

    lkp!(res, payload, lnRho, lnEi) = begin
        var, is_first = payload
        (is_first || isnothing(coefs[])) && (coefs[] = weights(eos, lnRho=lnRho, lnEi=lnEi))
        cRho, cT = coefs[]
        _dispatch_sample!(res, eos, var, cRho, cT)
    end

    (dispatch_eos_variables[q], q == first(dispatch_eos_quantities)), lkp!
end

dispatch_lookup_generator(eos::Union{RegularEoSTable,AxedEoS}, q) =
    dispatch_lookup_generator(extended(eos), q)

_dispatch_sample!(res, eos, var, cRho, cT) = begin
    sample!((res,), eos, (var,), cRho, cT)
    res .= exp.(res)
end
