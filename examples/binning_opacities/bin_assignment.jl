### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1fe66a04-11b2-11ee-1953-893c83adcc45
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using PythonPlot 
	using TSO
	using LaTeXStrings

	plt = pyplot;
end

# ╔═╡ 08a2cdc7-b0b6-4a37-bdce-b427452fd187
table = abspath("tables/TSO_MARCS_v1.6")

# ╔═╡ c7118fd1-a8b2-4672-88ed-bf3be87f51d8
binned_table = abspath("DIS_MARCS_v1.6.8") 

# ╔═╡ bd01db62-b3d7-4180-8995-53da826112a4
fopa = reload(
	SqOpacity, joinpath(table, "combined_formation_opacities.hdf5"), mmap=true
)

# ╔═╡ 4b1aeb85-ff6c-4436-a7ac-bc5907b15e66
assignment(path) = begin
	o = TSO.HDF5.h5open(path)
	TSO.HDF5.read(o["bins"])
end

# ╔═╡ 04685166-62c7-4944-b67c-09667329d546
a = assignment(joinpath(binned_table, "bin_assignment.hdf5"))

# ╔═╡ 9ef00b20-0707-4b34-9308-609ee898e409
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
		joinpath(binned_table, "formation_opacity.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═1fe66a04-11b2-11ee-1953-893c83adcc45
# ╠═08a2cdc7-b0b6-4a37-bdce-b427452fd187
# ╠═c7118fd1-a8b2-4672-88ed-bf3be87f51d8
# ╠═bd01db62-b3d7-4180-8995-53da826112a4
# ╠═4b1aeb85-ff6c-4436-a7ac-bc5907b15e66
# ╠═04685166-62c7-4944-b67c-09667329d546
# ╟─9ef00b20-0707-4b34-9308-609ee898e409
