using Pkg; Pkg.activate("..")
using TSO
using TimerOutputs

begin
    timer=TimerOutput()
    eos = reload("../../opacity_tables/magg_m0_a0_vmic1_v3.5/combined_eos_magg_m0_a0_vmic1.hdf5")
    opa = reload("../../opacity_tables/magg_m0_a0_vmic1_v3.5/combined_opacities_magg_m0_a0_vmic1.hdf5", mmap=true)
    
    lnt = range(3, 4.4, length=100) |> collect .|> exp10 .|> log
    lnrho = range(-12, -10, length=100) |> collect .|> exp10 .|> log
    
    reset_timer!(timer)
    lookup(eos, opa, :κ, lnrho, lnt)
    @timeit timer "old" k1 = lookup(eos, opa, :κ, lnrho, lnt)
    
    sample(eos, opa, [:κ], lnrho, lnt)
    @timeit timer "new" k2, = sample(eos, opa, [:κ], lnrho, lnt)
    show(timer)
    println("")
    @show eltype(k1) eltype(k2)
    @show typeof(k1) typeof(k2)
    @show size(k1) size(k2)
    d = abs.((permutedims(k1, (2, 1)) .- k2) ./ k2)
    @show maximum(d) minimum(d)

    k1, k2
end

