### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 80e0765c-5ecb-11ee-285f-0533135dd035
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot

	plt = pyplot; 
end

# ╔═╡ 04f9f6cd-b267-47b6-a901-6f48c6e88381
md"# Binning of Opacities"

# ╔═╡ 25f6b474-02e0-430e-98ea-dbf658b2199d
begin
	eos_folder = "../../../opacity_tables/TSO_M3D_magg_m0_a0_v1.4"
	model = joinpath(
		"../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5", "inim.dat"
	)
	name = "t5777g44m00"
	version = "v1.4"
    extension = "magg_m0_a0"
	eos_old_name = "combined_eos_$(extension).hdf5"
	opa_old_name = "combined_opacities_$(extension).hdf5"
end;

# ╔═╡ a0c3e42d-28aa-4b1e-8394-596833d53105
fopa_path = TSO.compute_formation_opacities(
	name, eos_folder, eos_old_name, opa_old_name, model, logg=4.40,
)
#fopa_path = joinpath(
#eos_folder, "combined_formation_opacities"*TSO.join_full(name, extension)*".hdf5"
#)

# ╔═╡ d19ec293-a019-4646-a316-82adbdd12a12
quadrantlimit(table_folder, name; extension, λ_lim=5.0) = begin
    name_ext = TSO.join_full(name, extension)
    ext = TSO.join_full(extension)

    opa = reload(
        SqOpacity, 
		joinpath(table_folder, "combined_opacities$(ext).hdf5"), 
		mmap=true
    )
	
    fopa = reload(
        SqOpacity, 
		joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"), 
		mmap=true
    )

    TSO.median(-log10.(fopa.κ_ross)[log10.(opa.λ) .> λ_lim])
end

# ╔═╡ a3c1263e-be42-4eb5-a1a0-c2213e9da771
ql = quadrantlimit(eos_folder, name, extension=extension, λ_lim=4.8)

# ╔═╡ 8691f075-0d53-4672-aa15-2cdf95e83980
new_eos_folder = TSO.bin_opacity_table(
	name,
	eos_folder,
	eos_old_name,
	opa_old_name,
	model,
	version=version,
	method=:kmeans, 
	Nbins=8, 
	quadrants=[ 
		TSO.Quadrant((0.0, 4.0), (ql, 4.5), 2, stripes=:κ),
		TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
		TSO.Quadrant((4.0, 100.0), (ql, 100), 1, stripes=:κ),
		TSO.Quadrant((0.0, 100.0), (-100, ql), 4, stripes=:λ),
	],
	maxiter=5000, display=:none
)

# ╔═╡ 40e3ea4c-c161-48ff-b9b6-a0736e229a5f
TSO.create_E_from_T(
	new_eos_folder, 
	upsample=2048, 
	version=version
)

# ╔═╡ 3c0c4f09-baad-4d33-aed2-b866c720861a
eos = reload(SqEoS, "DIS_MARCS_E_v1.0/eos.hdf5")

# ╔═╡ 10c779cb-445c-427d-96e2-ebe48ae93341
eos2 = reload(SqEoS, joinpath(eos_folder, eos_old_name))

# ╔═╡ 9b04ffb8-371d-457b-b5e3-cb360d97a4ca
md"# Visualization"

# ╔═╡ 2b79f6a7-3093-4c48-9e95-c6acdd3a2870
a = TSO.bin_assignment(joinpath(new_eos_folder, "bin_assignment.hdf5"))

# ╔═╡ d0f08fbd-503b-46dd-96b6-85b5c0879c7b
fopa = reload(SqOpacity, fopa_path)

# ╔═╡ abb97b43-be00-41fe-a7af-0d3ddabff04e
let
	plt.close()
	
	ff, axf = plt.subplots(figsize=(6,6))

	im = axf.scatter(
		log10.(fopa.λ), -log10.(fopa.κ_ross), s=0.7, cmap="rainbow"
	)	

	axf.set_ylabel(L"\rm \log\ \tau_{ross} (\tau_{\lambda}=1)", fontsize="large")
	axf.set_xlabel(L"\rm \lambda\ [\AA]", fontsize="large")
	
	axf.minorticks_on()
	axf.tick_params(top=true, right=true, direction="in", which="both")

	ff.colorbar(im, ax=axf)

	#ff.savefig(
	#	joinpath(new_eos_folder, "formation_opacity.png"), 
	#	dpi=800, bbox_inches="tight"
	#)
	
	gcf()
end

# ╔═╡ 39b1c56b-a7f8-4c05-8f72-5a792fa231c4
let
	f, ax = plt.subplots(2, 1, sharex=true)
	plt.subplots_adjust(hspace=0)

	m = TSO.flip!(Average3D(model))
	ax[0].plot(m.z, m.lnρ)
	ax[1].plot(m.z, m.lnT)

	ax[0].set_ylabel("lnρ")
	ax[1].set_ylabel("lnT")
	ax[1].set_xlabel("z")
	

	gcf()
end

# ╔═╡ 9c772bdc-b99d-4203-8b52-29752f48064a
let
	plt.close()

	ee, rr = TSO.meshgrid(@axed(eos))

	@show size(ee) size(rr) size(eos.lnEi)
	im = plt.scatter(ee, rr, c=eos.lnT)
	plt.colorbar(im)

	gcf()
end

# ╔═╡ 6459e4ae-7d7e-47e0-9e9c-08101573b36a
let
	plt.close()

	ee, rr = TSO.meshgrid(@axed(eos2))

	@show size(ee) size(rr) size(eos2.lnEi)
	im = plt.scatter(ee, log10.(exp.(rr)), c=exp.(eos2.lnEi))
	plt.colorbar(im)

	gcf()
end

# ╔═╡ Cell order:
# ╟─04f9f6cd-b267-47b6-a901-6f48c6e88381
# ╠═80e0765c-5ecb-11ee-285f-0533135dd035
# ╠═25f6b474-02e0-430e-98ea-dbf658b2199d
# ╠═a0c3e42d-28aa-4b1e-8394-596833d53105
# ╠═d19ec293-a019-4646-a316-82adbdd12a12
# ╠═a3c1263e-be42-4eb5-a1a0-c2213e9da771
# ╠═8691f075-0d53-4672-aa15-2cdf95e83980
# ╠═40e3ea4c-c161-48ff-b9b6-a0736e229a5f
# ╠═3c0c4f09-baad-4d33-aed2-b866c720861a
# ╠═10c779cb-445c-427d-96e2-ebe48ae93341
# ╟─9b04ffb8-371d-457b-b5e3-cb360d97a4ca
# ╠═2b79f6a7-3093-4c48-9e95-c6acdd3a2870
# ╠═d0f08fbd-503b-46dd-96b6-85b5c0879c7b
# ╠═abb97b43-be00-41fe-a7af-0d3ddabff04e
# ╟─39b1c56b-a7f8-4c05-8f72-5a792fa231c4
# ╠═9c772bdc-b99d-4203-8b52-29752f48064a
# ╠═6459e4ae-7d7e-47e0-9e9c-08101573b36a
