### A Pluto.jl notebook ###
# v0.19.30

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
	eos_folder = "tables/TSO_MARCS_v1.6"
	model = "sun_stagger.dat"
	name = "cnt"
	version = "v1.7.7"
    extension = "magg22"
end;

# ╔═╡ a0c3e42d-28aa-4b1e-8394-596833d53105
fopa_path = TSO.compute_formation_opacities(
	eos_folder, model, name, extension=extension, logg=4.438,
)
#fopa_path = joinpath(
#eos_folder, "combined_formation_opacities"*TSO.join_full(name, extension)*".hdf5"
#)

# ╔═╡ 8691f075-0d53-4672-aa15-2cdf95e83980
new_eos_folder = TSO.bin_opacity_table(
	eos_folder, 
	model, 
	name;
	version=version,
	extension=extension,
	method=:kmeans, 
	stripes=false,
	use_contribution=false, 
	Nbins=8, 
	quadrants=[ 
		TSO.Quadrant((0.0, 4.0), (1.0, 4.5), 2, stripes=:κ),
		TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
		TSO.Quadrant((4.0, 100.0), (1.0, 100), 1, stripes=:κ),
		TSO.Quadrant((0.0, 100.0), (-100, 1.0), 4, stripes=:λ),
	],
	maxiter=5000, display=:none
)

# ╔═╡ 40e3ea4c-c161-48ff-b9b6-a0736e229a5f
TSO.create_E_from_T(
	new_eos_folder, name, 
	upsample=2048, 
	version=version
)

# ╔═╡ 9b04ffb8-371d-457b-b5e3-cb360d97a4ca
md"# Visualization"

# ╔═╡ 2b79f6a7-3093-4c48-9e95-c6acdd3a2870
a = TSO.bin_assignment(joinpath(new_eos_folder, "bin_assignment.hdf5"))

# ╔═╡ d0f08fbd-503b-46dd-96b6-85b5c0879c7b
fopa = reload(SqOpacity, fopa_path)

# ╔═╡ abb97b43-be00-41fe-a7af-0d3ddabff04e
begin
	plt.close()
	
	ff, axf = plt.subplots(figsize=(6,6))

	im = axf.scatter(
		log10.(fopa.λ), -log10.(fopa.κ_ross), c=a, s=0.7, cmap="rainbow"
	)	

	axf.set_ylabel(L"\rm \log\ \tau_{ross} (\tau_{\lambda}=1)", fontsize="large")
	axf.set_xlabel(L"\rm \lambda\ [\AA]", fontsize="large")
	
	axf.minorticks_on()
	axf.tick_params(top=true, right=true, direction="in", which="both")

	ff.colorbar(im, ax=axf)

	ff.savefig(
		joinpath(new_eos_folder, "formation_opacity.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─04f9f6cd-b267-47b6-a901-6f48c6e88381
# ╠═80e0765c-5ecb-11ee-285f-0533135dd035
# ╠═25f6b474-02e0-430e-98ea-dbf658b2199d
# ╠═a0c3e42d-28aa-4b1e-8394-596833d53105
# ╠═8691f075-0d53-4672-aa15-2cdf95e83980
# ╠═40e3ea4c-c161-48ff-b9b6-a0736e229a5f
# ╟─9b04ffb8-371d-457b-b5e3-cb360d97a4ca
# ╠═2b79f6a7-3093-4c48-9e95-c6acdd3a2870
# ╠═d0f08fbd-503b-46dd-96b6-85b5c0879c7b
# ╟─abb97b43-be00-41fe-a7af-0d3ddabff04e
