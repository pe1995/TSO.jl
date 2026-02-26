using Pkg; Pkg.activate("..")
using TSO
using TimerOutputs

begin
    timer=TimerOutput()
    eos = reload("../../opacity_tables/magg_m0_a0_vmic1_v3.5/combined_eos_magg_m0_a0_vmic1.hdf5")
    opa = reload("../../opacity_tables/magg_m0_a0_vmic1_v3.5/combined_opacities_magg_m0_a0_vmic1.hdf5", mmap=true)
    
    lnt = range(3, 4.4, length=100) |> collect .|> exp10 .|> log
    lnrho = range(-12, -10, length=100) |> collect .|> exp10 .|> log
    
    #=reset_timer!(timer)
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
    @show maximum(d) minimum(d)=#

    # same test but with the new axis lookup
    reset_timer!(timer)
    lnpg_true, lnei_true = sample(eos, (:lnPg, :lnEi), lnrho, lnt)

    sample(eos, (:lnRho,), lnpg_true, lnt)
    sample(eos, (:lnT,), lnrho, lnei_true)
    sample(eos, opa, (:lnRho, :κ), lnpg_true, lnt)
    sample(eos, opa, (:lnT, :κ), lnrho, lnei_true)
    @timeit timer "lnrho_sampled_1" lnrho_sampled_1, = sample(eos, (:lnRho,), lnpg_true, lnt)
    @timeit timer "lnT_sampled_1" lnT_sampled_1, = sample(eos, (:lnT,), lnrho, lnei_true)
    #@timeit timer "lnrho_sampled_2" lnrho_sampled_2, _ = sample(eos, opa, (:lnRho, :κ), lnpg_true, lnt)
    #@timeit timer "lnT_sampled_2" lnT_sampled_2, _ = sample(eos, opa, (:lnT, :κ), lnrho, lnei_true)

    @show lnrho_sampled_1[99]
    @show lnrho[99]

    println("---------------")
    println("new:")
    
    d_rho = abs.((lnrho_sampled_1 .- lnrho) ./ lnrho)
    @show maximum(d_rho) minimum(d_rho)
    
    d_T = abs.((lnT_sampled_1 .- lnt) ./ lnt)
    @show maximum(d_T) minimum(d_T)

    #d_rho2 = abs.((lnrho_sampled_2 .- lnrho) ./ lnrho)
    #@show maximum(d_rho2) minimum(d_rho2)
    
    #d_T2 = abs.((lnT_sampled_2 .- lnt) ./ lnt)
    #@show maximum(d_T2) minimum(d_T2)

    #=println("---------------")
    println("old:")

    lnpg_true = lookup(eos, :lnPg, lnrho, lnt)
    lnei_true = lookup(eos, :lnEi, lnrho, lnt)
    #@timeit timer "old lnrho_sampled_1" lnrho_sampled_1 = lookup(eos, :lnRho, lnpg_true, lnt)
    lookup(eos, :lnT, lnrho, lnei_true)
    @timeit timer "old lnT_sampled_1" lnT_sampled_1 = lookup(eos, :lnT, lnrho, lnei_true)
    @timeit timer "old lnrho_sampled_1" lnrho_sampled_1 = [TSO.extended_lookup(extended(eos), :lnRho, lnpg_true[i], lnt[i]) for i in eachindex(lnpg_true)]


    
    d_rho = abs.((lnrho_sampled_1 .- lnrho) ./ lnrho)
    @show maximum(d_rho) minimum(d_rho)
    
    d_T = abs.((lnT_sampled_1 .- lnt) ./ lnt)
    @show maximum(d_T) minimum(d_T)=#

    println("")
    show(timer)
    println("")

    # test rosseland
    #=reset_timer!(timer)
    r1 = @timeit timer "old" TSO.rosseland_opacity_old(@axed(eos), opa)
    r2 = @timeit timer "new" rosseland_opacity(@axed(eos), opa)
    show(timer)
    println("")
    @show eltype(r1) eltype(r2)
    @show typeof(r1) typeof(r2)
    @show size(r1) size(r2)
    @show r1[100,100]
    @show r2[100,100]
    @show any(isnan, r1)
    @show any(isnan, r2)
    
    d = abs.((r1 .- r2) ./ r2)
    @show maximum(d) minimum(d)=#
end

