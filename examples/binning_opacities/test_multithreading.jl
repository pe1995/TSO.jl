begin
	using Pkg; Pkg.activate(".")
	using TSO
end

begin
	eos_folder = "../../../opacity_tables/TSO_M3D_magg_m0_a0_v2.0"
	av_path = joinpath(
		"../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5", "inim.dat"
	)
	name = "t5777g44m00"
	version = "t1.0"
    extension = "magg_m0_a0"
	eos_name = "combined_eos_$(extension).hdf5"
	opa_name = "combined_opacities_$(extension).hdf5"
	fopa_name = "combined_formation_opacities_$(name).hdf5"
end

begin
	# activate bin timing
	TSO.activate_timing!(TSO.binning_time)
	TSO.start_timing!()
end

begin
	eos = reload(SqEoS, joinpath(eos_folder, eos_name))
	opa = reload(SqOpacity, joinpath(eos_folder, opa_name), mmap=true)
	fopa = -log10.(reload(SqOpacity, joinpath(eos_folder, fopa_name), mmap=true).κ_ross)

	weights = ω_midpoint(opa)
	model = @optical TSO.flip(Average3D(eos, av_path, logg=4.44), depth=true) eos opa
	TSO.Clustering.Random.seed!(42)

	bin_assignment = ClusterBinning(
		TSO.KmeansBins;
		opacities=opa, 
		formation_opacity=fopa, 
		Nbins=8,
		quadrants=[ 
			TSO.Quadrant((0.0, 4.0), (1.0, 4.5), 2, stripes=:κ),
			TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
			TSO.Quadrant((4.0, 100.0), (1.0, 100), 1, stripes=:κ),
			TSO.Quadrant((0.0, 100.0), (-100, 1.0), 4, stripes=:λ),
		],
		maxiter=5000, 
		display=:none
	)

	bins = binning(bin_assignment, opa, fopa) 
end

begin
	binned_opacities = tabulate(bins, weights, eos, opa, transition_model=model)

	cp(joinpath(eos_folder, eos_name), "test_eos2.hdf5", force=true)
	save(binned_opacities.opacities, "test_binned_opacities2.hdf5")
end


begin
	# end bin timing
	TSO.end_timing!()
end