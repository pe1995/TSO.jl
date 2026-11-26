
function extend(tab::EquationOfState; downD=0.1, upD=0.0, downE=0.1, upE=0.0)
    ename = EnergyAxis(tab).name
    evals = EnergyAxis(tab).values

    dname = DensityAxis(tab).name
    dvals = DensityAxis(tab).values

    de = abs(maximum(evals) - minimum(evals))
    frac = abs((maximum(evals)+upE*de) - (minimum(evals)-downE*de)) / de
    evals_new = range(minimum(evals)-downE*de, maximum(evals)+upE*de, length=Int(ceil(frac * length(evals)))) |> collect

    de = abs(maximum(dvals) - minimum(dvals))
    frac = abs((maximum(dvals)+upD*de) - (minimum(dvals)-downD*de)) / de
    dvals_new = range(minimum(dvals)-downD*de, maximum(dvals)+upD*de, length=Int(ceil(frac * length(dvals)))) |> collect

    evals_new = Base.convert.(eltype(evals), evals_new)
    dvals_new = Base.convert.(eltype(dvals), dvals_new)
    interpolate_2D(tab; lnRho=dvals_new, (ename=>evals_new,)...)
end

