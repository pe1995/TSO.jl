#=
Extend EoS and opacity tables to values beyond the limits.
=#

"""
    extend(aos::AxedEoS, opacities...; downD=0.1, upD=0.0, downE=0.1, upE=0.0)

Extend the EoS upwards or downwards by `up` or `down` fractions of the current size.
The number of points will increase by the same fraction, however everything will
be interpolated to a new, uniform axis.
"""
function extend(aos::AxedEoS, opacities...; downD=0.1, upD=0.0, downE=0.1, upE=0.0)
    ename = EnergyAxis(aos).name
    evals = EnergyAxis(aos).values

    dname = DensityAxis(aos).name
    dvals = DensityAxis(aos).values

    # extend the axis up and down
    de = abs(maximum(evals) - minimum(evals))
    frac = abs((maximum(evals)+upE*de) - (minimum(evals)-downE*de)) / de
    evals_new = range(minimum(evals)-downE*de, maximum(evals)+upE*de, length=Int(ceil(frac * length(evals)))) |> collect

    de = abs(maximum(dvals) - minimum(dvals))
    frac = abs((maximum(dvals)+upD*de) - (minimum(dvals)-downD*de)) / de
    dvals_new = range(minimum(dvals)-downD*de, maximum(dvals)+upD*de, length=Int(ceil(frac * length(dvals)))) |> collect

    # convert to the correct type
    evals_new = Base.convert.(eltype(evals), evals_new)
    dvals_new = Base.convert.(eltype(dvals), dvals_new)
    interpolate_2D(aos, opacities...; lnRho=dvals_new, (ename=>evals_new,)...)
end

extend(eos::RegularEoSTable, args...; kwargs...) = extend(@axed(eos), args...; kwargs...)