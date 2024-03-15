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
	TSO.activate_timing!(TSO.rosseland_time)
	TSO.start_timing!()
end

begin
	eos = reload(SqEoS, joinpath(eos_folder, eos_name))
	opa = reload(SqOpacity, joinpath(eos_folder, opa_name), mmap=true)
end

begin
	eos_test = deepcopy(eos)
	TSO.rosseland_opacity!(eos_test.lnRoss, eos_test, opa)
end

begin
	save(eos_test, "eos_rosseland1.hdf5")
end

begin
	# end bin timing
	TSO.end_timing!()
end