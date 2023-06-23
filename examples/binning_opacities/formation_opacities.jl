### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 49b35982-11a8-11ee-3b21-2bdc91d7c4b5
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using PyPlot 
	using TSO
	using LaTeXStrings
end

# ╔═╡ 1f1efb64-f9e9-41ba-9796-2901ceda199d
table = abspath("tables/TSO_MARCS_v1.6")

# ╔═╡ fa9f67c9-7984-4593-8575-25cf9d370e53
table_old = abspath("tables/TSO_MARCS_v1.4")

# ╔═╡ a65d723f-1623-44cf-a8ec-90f53f32722c
fopa = reload(
	SqOpacity, joinpath(table, "combined_formation_opacities.hdf5"), mmap=true
)

# ╔═╡ a23d901d-e903-46d0-8106-c7a897839cab
begin
	ff, axf = plt.subplots(figsize=(5,6))

	axf.scatter(log10.(fopa.λ), -log10.(fopa.κ_ross), color="k", s=0.1)	

	axf.set_ylabel(L"\rm \log\ \tau_{ross} (\tau_{\lambda}=1)", fontsize="large")
	axf.set_xlabel(L"\rm \lambda\ [\AA]", fontsize="large")
	
	axf.minorticks_on()
	axf.tick_params(top=true, right=true, direction="in", which="both")


	ff.savefig(joinpath(table, "formation_opacity.png"), dpi=800, bbox_inches="tight")
	gcf()
end

# ╔═╡ Cell order:
# ╠═49b35982-11a8-11ee-3b21-2bdc91d7c4b5
# ╠═1f1efb64-f9e9-41ba-9796-2901ceda199d
# ╠═fa9f67c9-7984-4593-8575-25cf9d370e53
# ╠═a65d723f-1623-44cf-a8ec-90f53f32722c
# ╟─a23d901d-e903-46d0-8106-c7a897839cab
