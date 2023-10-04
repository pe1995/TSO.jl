### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 80e0765c-5ecb-11ee-285f-0533135dd035
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using Plots
end

# ╔═╡ 25f6b474-02e0-430e-98ea-dbf658b2199d
begin
	eos_folder = "tables/TSO_MARCS_v1.6"
	model = "sun_stagger.dat"
	name = "red"
	version = "v1.7.5"
    extension = "magg22"
end;

# ╔═╡ a0c3e42d-28aa-4b1e-8394-596833d53105
TSO.compute_formation_opacities(eos_folder, model, name, extension=extension)

# ╔═╡ 8691f075-0d53-4672-aa15-2cdf95e83980
new_eos_folder = TSO.bin_opacity_table(
	eos_folder, 
	model, 
	name;
	version=version,
	extension=extension,
	method=:kmeans, 
	use_contribution=false, 
	stripes=false,
	Nbins=7, 
	quadrants=[ 
		TSO.Quadrant((0.0, 4.0),   (-100, 4.5), 4),
		TSO.Quadrant((0.0, 4.0),   (4.5, 100), 4),
		TSO.Quadrant((4.0, 100.0), (-100, 100), 8)
	],
	maxiter=5000, display=:none
)

# ╔═╡ 40e3ea4c-c161-48ff-b9b6-a0736e229a5f
TSO.create_E_from_T(
	new_eos_folder, name, 
	upsample=2048, 
	version=version
)

# ╔═╡ Cell order:
# ╠═80e0765c-5ecb-11ee-285f-0533135dd035
# ╠═25f6b474-02e0-430e-98ea-dbf658b2199d
# ╠═a0c3e42d-28aa-4b1e-8394-596833d53105
# ╠═8691f075-0d53-4672-aa15-2cdf95e83980
# ╠═40e3ea4c-c161-48ff-b9b6-a0736e229a5f
