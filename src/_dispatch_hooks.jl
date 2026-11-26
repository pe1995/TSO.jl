# ============================================================================
# EoS hooks for the MUST.jl snapshot conversion
# ============================================================================

const dispatch_eos_quantities = [:T, :kr, :Pg, :Ne]
const dispatch_eos_variables = Dict(:T=>:lnT, :kr=>:lnRoss, :Pg=>:lnPg, :Ne=>:lnNe)

"""
    dispatch_eos_reader(eos_path; quantities=dispatch_eos_quantities)

`eos_reader` for `MUST.snapshotBox`, loading the internal energy grid of the
TSO table at `eos_path` (a table folder or file).
"""
dispatch_eos_reader(eos_path; quantities=dispatch_eos_quantities) = (run) -> begin
    eos = EquationOfState(eos_path; grid=:energy, parts=(:eos,))
    energy_variable(eos) == :lnEi || error(
        "$(eos_path) is not an energy based table, the conversion looks up in (lnRho, lnEi)."
    )

    eos, quantities
end

"""
    dispatch_lookup_generator(eos, q)

`lookup_generator` for `MUST.snapshotBox`. Returns the payload for quantity `q`
and the shared kernel that fills it in.
"""
function dispatch_lookup_generator(eos::EquationOfState, q)
    coefs = Ref{Any}(nothing)

    lkp!(res, payload, lnRho, lnEi) = begin
        var, is_first = payload
        (is_first || isnothing(coefs[])) && (coefs[] = weights(eos, lnRho=lnRho, lnEi=lnEi))
        cRho, cT = coefs[]
        _dispatch_sample!(res, eos, var, cRho, cT)
    end

    (dispatch_eos_variables[q], q == first(dispatch_eos_quantities)), lkp!
end

_dispatch_sample!(res, eos, var, cRho, cT) = begin
    sample!((res,), eos, (var,), cRho, cT)
    res .= exp.(res)
end
