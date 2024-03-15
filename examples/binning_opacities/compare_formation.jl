### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 3640be9e-e2d1-11ee-173c-e1e1b326e1af
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot
	plt = matplotlib.pyplot
end

# ╔═╡ 6836bbbf-a199-46b0-85ad-bf3f9da15172
fopa1 = reload(SqOpacity, "test_formation_opacities.hdf5")

# ╔═╡ ed9ca62a-54fd-44a0-a575-19f9740c5eef
fopa2 = reload(SqOpacity, "test_formation_opacities2.hdf5")

# ╔═╡ 357cf6ed-3152-4df0-b8cc-ec286ff22dda
let 
	plt.close()
	f, ax = plt.subplots(1, 2, sharey=true)

	ax[0].plot(log10.(fopa1.λ), -log10.(fopa1.κ_ross))
	ax[1].plot(log10.(fopa2.λ), -log10.(fopa2.κ_ross))
	
	gcf()
end

# ╔═╡ 0e5aacaa-6c4a-4035-a5e1-eb5b5f4caab8
let 
	plt.close()
	f, ax = plt.subplots(1, 1, sharey=true)

	ax.plot(log10.(fopa1.λ), -log10.(fopa1.κ_ross)+ log10.(fopa2.κ_ross))
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═3640be9e-e2d1-11ee-173c-e1e1b326e1af
# ╠═6836bbbf-a199-46b0-85ad-bf3f9da15172
# ╠═ed9ca62a-54fd-44a0-a575-19f9740c5eef
# ╠═357cf6ed-3152-4df0-b8cc-ec286ff22dda
# ╠═0e5aacaa-6c4a-4035-a5e1-eb5b5f4caab8
