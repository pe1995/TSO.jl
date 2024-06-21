### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ ce6b7716-2711-11ef-3dc4-9f6fb9c1c81b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using TSO
	using PythonPlot
end

# ╔═╡ 54e626e9-9fa3-4db8-bfe6-a9366c7d0c11
plt = matplotlib.pyplot

# ╔═╡ ca93fa86-b7b0-44b0-a1a2-ef0bbe5e50bd
matplotlib.style.use(joinpath(dirname(pathof(TSO)), "Bergemann2023.mplstyle"))

# ╔═╡ 93df09bf-c38e-4c62-9338-e077f6d48708


# ╔═╡ d99e4e3d-a867-4960-8f97-41936875d8bf
path = "/home/eitner/shared/StAt/opacity_tables/TSO_M3D_magg_m10_vmic2_v4.0"

# ╔═╡ af14a479-83f3-48fc-b15d-953ff4aa5286
eos = reload(SqEoS, joinpath(path, "combined_eos_magg_m10_vmic2.hdf5"))

# ╔═╡ 0ca38a6b-3545-4ec0-97a1-5d2a856710ef
opa = reload(SqOpacity, joinpath(path, "combined_opacities_magg_m10_vmic2.hdf5"), mmap=true)

# ╔═╡ 25e34b02-d70b-4e55-8f08-1820290b63f6


# ╔═╡ 059c507b-ce26-431d-9f93-f472de0646bf
path2 = "/home/eitner/shared/StAt/opacity_tables/TSO_MARCS_magg_m0_a0_v1.8/"

# ╔═╡ fab5b437-8e2a-4da8-83d9-ee74ac27a43f
eos2 = reload(SqEoS, joinpath(path2, "combined_eos_magg_m0_a0.hdf5"))

# ╔═╡ 9682cbd6-73c3-4b65-9483-d33a4c99e015
opa2 = reload(SqOpacity, joinpath(path2, "combined_opacities_magg_m0_a0.hdf5"), mmap=true)

# ╔═╡ 2098e07b-4cd5-414a-9134-d0aff7e9a6a8


# ╔═╡ 2c5bdaef-b59a-48b1-a33a-a58e9481e714
md"# Opacities"

# ╔═╡ a21b0082-8918-43e8-9e45-00ea61f16d04
T, lnρ = 6341.031935138665, -15.276813605421008

# ╔═╡ 05b762de-b559-4daf-a29e-4f3dd3e658c1
κ = log.(lookup(eos, opa, :κ, lnρ, log(T)))

# ╔═╡ 3768553b-8211-4c2a-b532-4125b2702aa1
κ2 = log.(lookup(eos2, opa2, :κ, lnρ, log(T)))

# ╔═╡ cfe93bed-843c-4085-bc16-82e1890d12cc
let
	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
	ax.plot(log10.(opa2.λ), κ2, lw=1, label="[Fe/H] = 0", color="k")
	ax.plot(log10.(opa.λ), κ, lw=2, label="[Fe/H] = -10", color="r")

	ax.set_xlabel("log10 λ [Å]")
	ax.set_ylabel("log10 κ [cm2 /g]")
	ax.legend()

	f
end

# ╔═╡ Cell order:
# ╠═ce6b7716-2711-11ef-3dc4-9f6fb9c1c81b
# ╠═54e626e9-9fa3-4db8-bfe6-a9366c7d0c11
# ╠═ca93fa86-b7b0-44b0-a1a2-ef0bbe5e50bd
# ╟─93df09bf-c38e-4c62-9338-e077f6d48708
# ╠═d99e4e3d-a867-4960-8f97-41936875d8bf
# ╠═af14a479-83f3-48fc-b15d-953ff4aa5286
# ╠═0ca38a6b-3545-4ec0-97a1-5d2a856710ef
# ╟─25e34b02-d70b-4e55-8f08-1820290b63f6
# ╠═059c507b-ce26-431d-9f93-f472de0646bf
# ╠═fab5b437-8e2a-4da8-83d9-ee74ac27a43f
# ╠═9682cbd6-73c3-4b65-9483-d33a4c99e015
# ╟─2098e07b-4cd5-414a-9134-d0aff7e9a6a8
# ╟─2c5bdaef-b59a-48b1-a33a-a58e9481e714
# ╠═a21b0082-8918-43e8-9e45-00ea61f16d04
# ╠═05b762de-b559-4daf-a29e-4f3dd3e658c1
# ╠═3768553b-8211-4c2a-b532-4125b2702aa1
# ╟─cfe93bed-843c-4085-bc16-82e1890d12cc
