### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ c396ea3c-e215-11ee-3d39-af55ebc8b47c
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using PythonPlot
	using TSO
	plt = matplotlib.pyplot
end

# ╔═╡ e69c945c-280f-4644-9daa-26b95bf80288
eos1 = reload(SqEoS, "test_eos.hdf5")

# ╔═╡ 6e323b54-fce1-4eab-bb0f-fed855a56ad5
opa1 = reload(SqOpacity, "test_binned_opacities.hdf5")

# ╔═╡ 1933754c-d118-4af1-8d95-d5c756ec7b92
eos2 = reload(SqEoS, "test_eos2.hdf5")

# ╔═╡ 43904148-3a76-46a6-b184-11ebbb0e80c6
opa2 = reload(SqOpacity, "test_binned_opacities2.hdf5")

# ╔═╡ 6d5045f9-9052-478a-8908-81cdc5c778cd
rosseland2 = reload(SqEoS, "eos_rosseland.hdf5")

# ╔═╡ 55ae5dd1-2bca-4a00-93f2-971a0dac9a9b
rosseland1 = reload(SqEoS, "eos_rosseland1.hdf5")

# ╔═╡ 0629c0d2-aa1f-4e81-8532-a191756fbdaa
T, ρ = 6000., 1e-6

# ╔═╡ 1264403c-746e-4def-8719-918aea46bab4
let
	plt.close()

	k1 = log10.(lookup(eos1, opa1, :κ, log(ρ), log(T)))
	k2 = log10.(lookup(eos2, opa2, :κ, log(ρ), log(T)))
	
	plt.plot(opa1.λ, k1, label="1 thread", lw=4, color="black")
	plt.plot(opa2.λ, k2, label="multi thread", lw=2, color="magenta", ls="--")

	plt.ylabel("log opacity")
	plt.xlabel("bin")

	plt.legend()
	
	gcf()
end

# ╔═╡ fd70d88d-a3d1-4f8d-832a-fdbe651ceb9f
let
	plt.close()

	f, ax = plt.subplots(1, 2)
	xx, yy = meshgrid(@axed eos1)
	ax[0].scatter(xx, yy, c=rosseland1.lnRoss)
	ax[1].scatter(xx, yy, c=rosseland2.lnRoss)

	gcf()
end

# ╔═╡ Cell order:
# ╠═c396ea3c-e215-11ee-3d39-af55ebc8b47c
# ╠═e69c945c-280f-4644-9daa-26b95bf80288
# ╠═6e323b54-fce1-4eab-bb0f-fed855a56ad5
# ╠═1933754c-d118-4af1-8d95-d5c756ec7b92
# ╠═43904148-3a76-46a6-b184-11ebbb0e80c6
# ╠═6d5045f9-9052-478a-8908-81cdc5c778cd
# ╠═55ae5dd1-2bca-4a00-93f2-971a0dac9a9b
# ╠═0629c0d2-aa1f-4e81-8532-a191756fbdaa
# ╟─1264403c-746e-4def-8719-918aea46bab4
# ╠═fd70d88d-a3d1-4f8d-832a-fdbe651ceb9f
