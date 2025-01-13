### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 22d20ce6-e5f1-11ee-0954-110219aa72b9
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot
	using LaTeXStrings
end

# ╔═╡ 414bbbfb-4eb5-4767-9509-f750c168e2c3
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(TSO)), "Bergemann2023.mplstyle"))
end

# ╔═╡ 17204fec-a78a-424d-ad0f-c587eb063f64


# ╔═╡ 29c10911-c084-4811-8cec-e2466ac1826d
md"# Unbinned Opacities"

# ╔═╡ 72f54a34-8bfc-4dd8-a3bf-0932f4451e41
#eosUnbinnedPath = "../../../opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"
eosUnbinnedPath = "../../../opacity_tables/magg_m0_a0_vmic1_v1.1"

# ╔═╡ 837c3da9-7618-481e-b692-c4bd8217c429
#=eosUnbinned = reload(
	SqEoS, 
	joinpath(eosUnbinnedPath, "combined_eos_magg_m0_a0.hdf5")
)=#
eosUnbinned = reload(
	SqEoS, 
	joinpath(eosUnbinnedPath, "combined_eos_magg_m0_a0_vmic1.hdf5")
)

# ╔═╡ 235df4e6-088b-48a2-8658-78a5b894d2d0
begin
	#=opaUnbinned = reload(
		SqOpacity, 
		joinpath(eosUnbinnedPath, "combined_opacities_magg_m0_a0.hdf5"), 
		mmap=true
	)=#
	opaUnbinned = reload(
		SqOpacity, 
		joinpath(eosUnbinnedPath, "combined_opacities_magg_m0_a0_vmic1.hdf5"), 
		mmap=false
	)
	T = exp.(eosUnbinned.lnT)
	for j in eachindex(eosUnbinned.lnRho)
		opaUnbinned.src[:, j, :] .= TSO.Bλ(opaUnbinned.λ, T)
	end
end

# ╔═╡ 42e48cc3-77e6-475e-bdf3-b942fd28106f


# ╔═╡ 16af61db-c470-4634-93b5-67d94ef6983d


# ╔═╡ c8fd5b7e-6dce-44a5-a8a6-6c1c84dba0d2
md"# Binned Opacities"

# ╔═╡ e9f47028-8d85-4d46-b12b-d1bff9489275
eosBinnedPath = "../../../stellar_atmospheres/MainSequenceInterpolated/ST22D6_E_t46.26g34.98m0.000_v1.0"
#eosBinnedPath = "../../../stellar_atmospheres/MainSequenceInterpolated/ST22C_E_t57.77g44.40m0.000_v1.0"

# ╔═╡ 73a19514-66f8-4508-b561-413592380098
eosBinned = reload(
	SqEoS, 
	joinpath(eosBinnedPath, "eos_T.hdf5")
)

# ╔═╡ 6244bbdd-bf0c-4d2f-bf8d-ba5b6c0daf1a
opaBinned = reload(
	SqOpacity, 
	joinpath(eosBinnedPath, "binned_opacities_T.hdf5")
)

# ╔═╡ 002d9db6-c483-4422-957c-86e1802b551f


# ╔═╡ 4b20e231-420b-4813-ae9e-8f439773ed78


# ╔═╡ 23fa569d-6fb4-4dd8-9a62-e7f2abd85de3
md"# Model Atmosphere"

# ╔═╡ 9ef995ba-cd4a-4ab0-98d8-42a68f832dd6
model = TSO.@optical(
	TSO.Average3D(eosUnbinned, joinpath(eosBinnedPath, "inim.dat")), eosUnbinned, opaUnbinned
)
#=model = TSO.@optical(
	TSO.Average3D(
		eosUnbinned, 
		"../../../MUST.jl/examples/initial_models/av_models/t65g45m00_00069_av"
	), 
	eosUnbinned, 
	opaUnbinned
)=#

# ╔═╡ 83619cb0-6edb-46ff-8367-0378b2b2783e


# ╔═╡ da26dc03-a2c2-4333-a0ee-403039d20c65
md"# Radiative Transfer"

# ╔═╡ 4dbade62-2d8d-46b0-83ad-1c3982ddf5c9
weights = TSO.ω_midpoint(opaUnbinned)

# ╔═╡ 7db50062-1c29-4768-a9f5-7a24259850e3


# ╔═╡ 7e657d06-8fa1-4fc2-8427-3293701bff60
md"For the transfer solver we have to wrapped the binned opacities, so that the solver recognizes them as binned, because opacity is multiplied by ρ."

# ╔═╡ 2049e27b-bbc0-4e33-9e2c-f1f0842e97f3
opacityBinned = TSO.@binned opaBinned

# ╔═╡ 05b6c979-f636-4efe-9497-80237b9f2001


# ╔═╡ bf9fd12b-6361-4188-8bf4-1db0fc001679
md"Now the same can be done for the unbinned opacities."

# ╔═╡ db154864-5c13-46f4-8460-3764ec871f40
opacityUnbinned = TSO.@binned opaUnbinned eosUnbinned

# ╔═╡ 49a21c62-8ff0-4a80-93e0-b82f801f6d00


# ╔═╡ 513808b0-bb41-40b4-b805-a2e33dcaf8b5
md"From those we can now construct the solvers and compute the heating rates. Make sure to include the frequency weights in the unbinned case, which are already included in the binned opacities."

# ╔═╡ 11df28d5-73a7-4bd0-9303-baf84f916668
solverBinned = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinned), 
	opacities=opacityBinned
)

# ╔═╡ 4b4da494-897c-4a21-92eb-65f59901fec0
solverUnbinned = TSO.Heating(
	model, 
	eos=TSO.@axed(eosUnbinned), 
	opacities=opacityUnbinned
)

# ╔═╡ 37c19400-6f85-4141-ba1b-d70e7679a2be


# ╔═╡ 32ed8ad6-a66f-4e70-ae24-76d5f1a34224


# ╔═╡ 60fd07be-18fe-43bb-934a-8467a6753fd5
QBinned = TSO.solve(solverBinned)

# ╔═╡ 284edea4-bb09-4efb-9f5a-de9bda495914
QUnbinned = TSO.solve(solverUnbinned, weights=weights)

# ╔═╡ db1b70dc-9f0c-473d-94fc-2d87dacc1bc9


# ╔═╡ 1b039417-c067-4ec6-9047-7b1f27ca16ef
md"# Comparison"

# ╔═╡ 58498690-62ba-46e9-955c-f964970cb052
z, lnT, τ = QUnbinned.model.z, QUnbinned.model.lnT, QUnbinned.model.τ

# ╔═╡ 60138269-4108-44a4-9b4d-dcf397931be0
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	norm = maximum(abs.(TSO.heating(QUnbinned)))
	ax.plot(
		log10.(τ), TSO.heating(QUnbinned)./norm, 
		marker="x", color="k", ls="",
		markerfacecolor="w", markersize=7, label="unbinned"
	)
	
	ax.plot(
		log10.(τ), TSO.heating(QBinned)./norm, 
		marker="", color="k", ls="-",
		label="binned", lw=1
	)

	ax.set_ylabel(L"\rm q_{rad}\ /\ \left|q_{rad}^{\lambda, max}\right|")	
	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_xlim(-6, 4)

	ax.legend()
	
	gcf()
end

# ╔═╡ 56acf259-9998-4d39-a99b-2bf4cfaf37a6
let
	plt.close()

	f, ax = plt.subplots(figsize=(5,6))

	norm = maximum(abs.(TSO.heating(QUnbinned)))
	mask = sortperm(TSO.heating(QUnbinned))

	diag = range(
		-1.1,
		0.1,
		length=1000
	)
	
	ax.plot(diag, diag, color="k", ls=":")

	ax.set_xlabel(L"\rm q_{rad}^{\lambda}\ /\  \left|q_{rad}^{\lambda, max}\right|")
	ax.set_ylabel(L"\rm q_{rad}^{bin}\ /\  \left|q_{rad}^{\lambda, max}\right|")
	ax.set_xlim(-1.1, 0.1)
	ax.set_ylim(-1.1, 0.1)
	
	ax.legend()
	
	gcf()
end

# ╔═╡ f461a848-0d30-4fb6-81ca-27905544d4e0
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	norm = maximum(abs.(TSO.heating(QUnbinned))) / 100.
	stat_label(yi) = "$(round(maximum(abs.(y))/norm, sigdigits=3)) %"

	
	
	
	y = (TSO.heating(QBinned) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="cyan", label=L"\rm max(\delta q) = "*stat_label(y),
		lw=3, ls="-"
	)

	
	ax.set_ylabel(L"\rm \left(q_{rad}^{bin} - q_{rad}^{\lambda} \right)\ /\ \left|q_{rad}^{\lambda, max}\right|\ [\%]")	
	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_xlim(-4, 4)
	ax.set_ylim(-3.5, 3.5)

	ax.legend(
		labelspacing=0.05, handlelength=2, loc="lower center", ncols=2,
		bbox_to_anchor=(0.5, 0.99)
	)
	
	gcf()
end

# ╔═╡ 5b5243b0-189f-48ca-a8c2-5aa63f7c74c1


# ╔═╡ Cell order:
# ╠═22d20ce6-e5f1-11ee-0954-110219aa72b9
# ╠═414bbbfb-4eb5-4767-9509-f750c168e2c3
# ╟─17204fec-a78a-424d-ad0f-c587eb063f64
# ╟─29c10911-c084-4811-8cec-e2466ac1826d
# ╠═72f54a34-8bfc-4dd8-a3bf-0932f4451e41
# ╠═837c3da9-7618-481e-b692-c4bd8217c429
# ╠═235df4e6-088b-48a2-8658-78a5b894d2d0
# ╟─42e48cc3-77e6-475e-bdf3-b942fd28106f
# ╟─16af61db-c470-4634-93b5-67d94ef6983d
# ╟─c8fd5b7e-6dce-44a5-a8a6-6c1c84dba0d2
# ╠═e9f47028-8d85-4d46-b12b-d1bff9489275
# ╟─73a19514-66f8-4508-b561-413592380098
# ╟─6244bbdd-bf0c-4d2f-bf8d-ba5b6c0daf1a
# ╟─002d9db6-c483-4422-957c-86e1802b551f
# ╟─4b20e231-420b-4813-ae9e-8f439773ed78
# ╟─23fa569d-6fb4-4dd8-9a62-e7f2abd85de3
# ╠═9ef995ba-cd4a-4ab0-98d8-42a68f832dd6
# ╟─83619cb0-6edb-46ff-8367-0378b2b2783e
# ╟─da26dc03-a2c2-4333-a0ee-403039d20c65
# ╠═4dbade62-2d8d-46b0-83ad-1c3982ddf5c9
# ╟─7db50062-1c29-4768-a9f5-7a24259850e3
# ╟─7e657d06-8fa1-4fc2-8427-3293701bff60
# ╟─2049e27b-bbc0-4e33-9e2c-f1f0842e97f3
# ╟─05b6c979-f636-4efe-9497-80237b9f2001
# ╟─bf9fd12b-6361-4188-8bf4-1db0fc001679
# ╟─db154864-5c13-46f4-8460-3764ec871f40
# ╟─49a21c62-8ff0-4a80-93e0-b82f801f6d00
# ╟─513808b0-bb41-40b4-b805-a2e33dcaf8b5
# ╟─11df28d5-73a7-4bd0-9303-baf84f916668
# ╟─4b4da494-897c-4a21-92eb-65f59901fec0
# ╟─37c19400-6f85-4141-ba1b-d70e7679a2be
# ╟─32ed8ad6-a66f-4e70-ae24-76d5f1a34224
# ╟─60fd07be-18fe-43bb-934a-8467a6753fd5
# ╠═284edea4-bb09-4efb-9f5a-de9bda495914
# ╟─db1b70dc-9f0c-473d-94fc-2d87dacc1bc9
# ╟─1b039417-c067-4ec6-9047-7b1f27ca16ef
# ╠═58498690-62ba-46e9-955c-f964970cb052
# ╟─60138269-4108-44a4-9b4d-dcf397931be0
# ╟─56acf259-9998-4d39-a99b-2bf4cfaf37a6
# ╟─f461a848-0d30-4fb6-81ca-27905544d4e0
# ╟─5b5243b0-189f-48ca-a8c2-5aa63f7c74c1
