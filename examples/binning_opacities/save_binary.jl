### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 7acbc79a-222c-11ee-3455-79375786e57b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
end

# ╔═╡ 7e754ba0-1e0b-44c5-9cf0-9e564ca6c9b5
eos = reload(SqEoS, "binned/DIS_MARCS_v1.7.2/eos.hdf5")

# ╔═╡ 2d034698-70ca-40eb-95dc-e82b7e304ba9
opa = reload(SqOpacity, "binned/DIS_MARCS_v1.7.2/binned_opacities.hdf5")

# ╔═╡ bd63928c-9b53-43c7-be44-5c007c464b69
open("magg_12.txt", "w") do f
	write(f, "$(length(eos.lnT))\n")
	write(f, "$(length(eos.lnRho))\n")
	write(f, "$(length(opa.λ))")
end

# ╔═╡ c8344d0d-d79e-454e-ac33-dfca3d369ea8
open("magg_12.bin", "w") do f
	write(f, eos.lnT)
	write(f, eos.lnRho)
	write(f, opa.κ)
	write(f, opa.src)
end

# ╔═╡ 58d3e9f5-09c2-4d31-92e9-7e2ec3d2f6bf
eos.lnT

# ╔═╡ Cell order:
# ╠═7acbc79a-222c-11ee-3455-79375786e57b
# ╠═7e754ba0-1e0b-44c5-9cf0-9e564ca6c9b5
# ╠═2d034698-70ca-40eb-95dc-e82b7e304ba9
# ╠═bd63928c-9b53-43c7-be44-5c007c464b69
# ╠═c8344d0d-d79e-454e-ac33-dfca3d369ea8
# ╠═58d3e9f5-09c2-4d31-92e9-7e2ec3d2f6bf
