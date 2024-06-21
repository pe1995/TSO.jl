### A Pluto.jl notebook ###
# v0.19.38

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
eosUnbinnedPath = "../../../opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"

# ╔═╡ 837c3da9-7618-481e-b692-c4bd8217c429
eosUnbinned = reload(
	SqEoS, 
	joinpath(eosUnbinnedPath, "ross_combined_eos_magg_m0_a0.hdf5")
)

# ╔═╡ 235df4e6-088b-48a2-8658-78a5b894d2d0
opaUnbinned = reload(
	SqOpacity, 
	joinpath(eosUnbinnedPath, "combined_opacities_magg_m0_a0.hdf5"), 
	mmap=true
)

# ╔═╡ 16af61db-c470-4634-93b5-67d94ef6983d


# ╔═╡ c8fd5b7e-6dce-44a5-a8a6-6c1c84dba0d2
md"# Binned Opacities"

# ╔═╡ e9f47028-8d85-4d46-b12b-d1bff9489275
eosBinnedPath = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5"

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

# ╔═╡ 5aa6c145-f149-49ae-8644-64d5d6ae7ac1


# ╔═╡ 506d9c37-3e08-477c-a6b8-dc735c271672
eosBinnedMuramPath = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5.2"

# ╔═╡ 9ee859b8-84b4-4ff9-84a9-6c20244313bb
eosBinnedMuram = reload(
	SqEoS, 
	joinpath(eosBinnedMuramPath, "eos_T.hdf5")
)

# ╔═╡ 1ec7f614-aabd-46d2-b392-c3c3f2c5093a
opaBinnedMuram = reload(
	SqOpacity, 
	joinpath(eosBinnedMuramPath, "binned_opacities_T.hdf5")
)

# ╔═╡ 0e51622e-0f0b-41db-acee-487ea486d75f


# ╔═╡ 1b31714e-6b07-4cee-9205-a907a4fae9ef
eosBinnedPaperPath = "../../../opacity_tables/DIS_MARCS_v1.6.3"

# ╔═╡ 4d08803b-9e83-42f7-aec9-f2b4d0fca661
eosBinnedPaper = reload(
	SqEoS, 
	joinpath(eosBinnedPaperPath, "eos.hdf5")
)

# ╔═╡ e2b59fe6-0a6a-4980-b390-79699db977ca
opaBinnedPaper = reload(
	SqOpacity, 
	joinpath(eosBinnedPaperPath, "binned_opacities.hdf5")
)

# ╔═╡ 2eb2a0c0-fdae-4290-994b-8dc39b9f88ca


# ╔═╡ eb907c55-6357-49b9-b022-542a24741585
eosBinnedPaper5Path = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5.3"

# ╔═╡ 2211cd3c-1a8b-49b7-ad39-1c7849a5875e
eosBinnedPaper5 = reload(
	SqEoS, 
	joinpath(eosBinnedPaper5Path, "eos.hdf5")
)

# ╔═╡ eb4669a5-aecb-4320-8477-9ef581a6db2f
opaBinnedPaper5 = reload(
	SqOpacity, 
	joinpath(eosBinnedPaper5Path, "binned_opacities.hdf5")
)

# ╔═╡ 1f568201-92d4-493e-a9d0-7d8ca42caeb9
eosBinnedPaper7Path = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5.6"

# ╔═╡ 8ce20d39-2a0c-490c-97f9-2e21559a6bdf
eosBinnedPaper7 = reload(
	SqEoS, 
	joinpath(eosBinnedPaper7Path, "eos.hdf5")
)

# ╔═╡ 3dd054b6-65a8-48c5-a91c-14010aa06725
opaBinnedPaper7 = reload(
	SqOpacity, 
	joinpath(eosBinnedPaper7Path, "binned_opacities.hdf5")
)

# ╔═╡ ff459b9c-41c4-4ef4-802c-fd486b5e1504


# ╔═╡ 24dc0315-e841-4eb6-8d3b-b8c6ba5c91f8
eosBinnedPaper8Path = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5.4"

# ╔═╡ 76dcf498-0fe5-467c-bbb7-ed66d34e7f46
eosBinnedPaper8 = reload(
	SqEoS, 
	joinpath(eosBinnedPaper8Path, "eos.hdf5")
)

# ╔═╡ a833b314-b7a3-4bf4-a092-deb9875d3a7b
opaBinnedPaper8 = reload(
	SqOpacity, 
	joinpath(eosBinnedPaper8Path, "binned_opacities.hdf5")
)

# ╔═╡ f92403c6-75bf-4886-bd96-96cd8aa94146


# ╔═╡ f3650b49-a046-4caf-b97f-51a476764fff
eosBinnedPaper12Path = "../../../MUST.jl/examples/initial_models/DIS_MARCS_E_t5777g44m00_v0.5.5"

# ╔═╡ 02130696-3615-42c2-a4c0-5d15a54d2f3a
eosBinnedPaper12 = reload(
	SqEoS, 
	joinpath(eosBinnedPaper12Path, "eos.hdf5")
)

# ╔═╡ b824b498-c3ec-481a-89d7-922bcc9e273f
opaBinnedPaper12 = reload(
	SqOpacity, 
	joinpath(eosBinnedPaper12Path, "binned_opacities.hdf5")
)

# ╔═╡ 4b20e231-420b-4813-ae9e-8f439773ed78


# ╔═╡ 23fa569d-6fb4-4dd8-9a62-e7f2abd85de3
md"# Model Atmosphere"

# ╔═╡ 9ef995ba-cd4a-4ab0-98d8-42a68f832dd6
model = TSO.@optical(
	TSO.Average3D(eosUnbinned, "sun_stagger.dat"), eosUnbinned, opaUnbinned
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

# ╔═╡ 0c64314d-3c2b-43df-85e1-121e3ad2c2ba
opacityBinnedMuram = TSO.@binned opaBinnedMuram

# ╔═╡ a6d9cb7c-28f1-48ad-8ff3-aaa872ef5da8
opacityBinnedPaper = TSO.@binned opaBinnedPaper

# ╔═╡ 611fbf7c-51de-457f-a35f-230db8992979
opacityBinnedPaper5 = TSO.@binned opaBinnedPaper5

# ╔═╡ d1848d82-a027-4e1d-bdff-9a60a8a1986b
opacityBinnedPaper7 = TSO.@binned opaBinnedPaper7

# ╔═╡ 9f9ab022-e122-4880-b2a5-419def02de86
opacityBinnedPaper8 = TSO.@binned opaBinnedPaper8

# ╔═╡ 7dc65fe5-575f-458f-a562-26f0adc380c0
opacityBinnedPaper12 = TSO.@binned opaBinnedPaper12

# ╔═╡ aed8567f-b982-4da6-a913-69091fe49962


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

# ╔═╡ 43a0bf36-ed3d-49af-9be2-a26295a977f6
solverBinnedMuram = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedMuram), 
	opacities=opacityBinnedMuram
)

# ╔═╡ 6e9f815b-fd56-4bd3-a2d0-ccc89ee81c64
solverBinnedPaper = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedPaper), 
	opacities=opacityBinnedPaper
)

# ╔═╡ d9b1484a-c91c-4def-956c-3fb05469c6ce
solverBinnedPaper5 = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedPaper5), 
	opacities=opacityBinnedPaper5
)

# ╔═╡ d524e5fb-0622-4503-8951-dbfc2be84088
solverBinnedPaper7 = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedPaper7), 
	opacities=opacityBinnedPaper7
)

# ╔═╡ 8ba7a253-482e-44f8-b7e3-5a3485bf343e
solverBinnedPaper8 = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedPaper8), 
	opacities=opacityBinnedPaper8
)

# ╔═╡ f2497cb5-efd6-4045-80d2-63691055a70b
solverBinnedPaper12 = TSO.Heating(
	model, 
	eos=TSO.@axed(eosBinnedPaper12), 
	opacities=opacityBinnedPaper12
)

# ╔═╡ e6d8b403-f8d0-444e-a994-e927d54b7d80


# ╔═╡ 4b4da494-897c-4a21-92eb-65f59901fec0
solverUnbinned = TSO.Heating(
	model, 
	eos=TSO.@axed(eosUnbinned), 
	opacities=opacityUnbinned
)

# ╔═╡ 37c19400-6f85-4141-ba1b-d70e7679a2be


# ╔═╡ 60fd07be-18fe-43bb-934a-8467a6753fd5
QBinned = TSO.solve(solverBinned)

# ╔═╡ 54bbc506-0c18-4ab3-a3c7-1d811c54f73e
QBinnedMuram = TSO.solve(solverBinnedMuram)

# ╔═╡ 5a600020-6012-4ba4-8ed2-cd08ec7bc4d5
QBinnedPaper = TSO.solve(solverBinnedPaper)

# ╔═╡ 13f694d0-9ea8-413c-943a-b196baad5d89
QBinnedPaper5 = TSO.solve(solverBinnedPaper5)

# ╔═╡ 39b151a7-7d0b-4474-9c94-289136fe4c50
QBinnedPaper7 = TSO.solve(solverBinnedPaper7)

# ╔═╡ 07644fb0-3db1-4e46-adf6-5bac7cdb3bf3
QBinnedPaper8 = TSO.solve(solverBinnedPaper8)

# ╔═╡ d7951d3a-2d7b-4303-aa3d-0b326ea51e53
QBinnedPaper12 = TSO.solve(solverBinnedPaper12)

# ╔═╡ 38fe75d7-ee09-4dbc-aedd-fbf8ef304bb9


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
		log10.(τ), TSO.heating(QBinnedPaper5)./norm, 
		color="cyan", label="5 bins", lw=3
	)
	ax.plot(
		log10.(τ), TSO.heating(QBinnedPaper)./norm, 
		color="k", label="7 bins", lw=3
	)
	ax.plot(
		log10.(τ), TSO.heating(QBinnedPaper8)./norm, 
		color="magenta", label="8 bins", lw=3
	)
	ax.plot(
		log10.(τ), TSO.heating(QBinnedPaper12)./norm, 
		color="lime", label="12 bins", lw=2, ls="--"
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
	ax.plot(
		TSO.heating(QUnbinned)[mask]./norm, TSO.heating(QBinnedPaper5)[mask]./norm, 
		color="cyan", label="5 bins", lw=1, ls="-", marker="s", markersize=9
	)
	ax.plot(
		TSO.heating(QUnbinned)[mask]./norm, TSO.heating(QBinnedPaper)[mask]./norm, 
		color="k", label="7 bins", lw=1, marker="v", ls="-", markersize=9
	)
	ax.plot(
		TSO.heating(QUnbinned)[mask]./norm, TSO.heating(QBinnedPaper8)[mask]./norm, 
		color="magenta", label="8 bins", lw=1, marker=".", ls="-", markersize=9
	)
	ax.plot(
		TSO.heating(QUnbinned)[mask]./norm, TSO.heating(QBinnedPaper12)[mask]./norm, 
		color="lime", label="12 bins", lw=1, marker="x", ls="-", markersize=9
	)

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

	
	
	y = (TSO.heating(QBinnedPaper5) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="k", label="5 bins, "*L"\rm max(\delta q) = "*stat_label(y), 
		lw=2, ls=":"
	)
	
	y = (TSO.heating(QBinnedPaper7) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="magenta", label="7 bins, "*L"\rm max(\delta q) = "*stat_label(y),
		lw=3, ls="-"
	)

	y = (TSO.heating(QBinnedPaper8) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ),y./norm, 
		color="k", label="8 bins, "*L"\rm max(\delta q) = "*stat_label(y),
		lw=2, ls="--"
	)
	
	y = (TSO.heating(QBinnedPaper12) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="k", label="12 bins, "*L"\rm max(\delta q) = "*stat_label(y),
		lw=2, ls="-."
	)

	y = (TSO.heating(QBinned) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="cyan", label="8 bins (new method), "*L"\rm max(\delta q) = "*stat_label(y),
		lw=3, ls="-"
	)

	y = (TSO.heating(QBinnedMuram) .- TSO.heating(QUnbinned))
	ax.plot(
		log10.(τ), y./norm, 
		color="lime", label="MURaM bins, "*L"\rm max(\delta q) = "*stat_label(y), 
		lw=2, ls="-"
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


# ╔═╡ 87a57724-8403-44ed-85ef-6d101409b299
md"# Assignment"

# ╔═╡ 61b38b1c-e9be-43f6-a642-be6b72a27d5c
bin_assignment(path) = begin
    fid = TSO.HDF5.h5open(path, "r")
    bins = TSO.HDF5.read(fid["bins"])
    λ = TSO.HDF5.read(fid["lambda"])
    close(fid)

    bins, λ
end

# ╔═╡ 50e0e483-82f9-4b58-baac-ea0df32b8cfa
a, lam = bin_assignment(joinpath(eosBinnedPath, "bin_assignment.hdf5"))

# ╔═╡ 5f9d2291-2d6b-4a0a-a19d-c3a9ac28707e
fopa = reload(
	SqOpacity, 
	joinpath(eosUnbinnedPath, "combined_formation_opacities_t5777g44m00.hdf5"), 
	mmap=true
)

# ╔═╡ de186f85-91b3-4dd8-9f66-aca0ade1cb03
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	im = ax.scatter(
		log10.(fopa.λ), -log10.(fopa.κ_ross),
		s=1, c=a, cmap="rainbow"
	)
	f.colorbar(im, ax=ax)

	ax.set_ylabel(L"\rm -\log \tau_{ross}\left(\ \tau_{\lambda}=1\ \right)")

	gcf()
end

# ╔═╡ Cell order:
# ╠═22d20ce6-e5f1-11ee-0954-110219aa72b9
# ╠═414bbbfb-4eb5-4767-9509-f750c168e2c3
# ╟─17204fec-a78a-424d-ad0f-c587eb063f64
# ╟─29c10911-c084-4811-8cec-e2466ac1826d
# ╠═72f54a34-8bfc-4dd8-a3bf-0932f4451e41
# ╟─837c3da9-7618-481e-b692-c4bd8217c429
# ╟─235df4e6-088b-48a2-8658-78a5b894d2d0
# ╟─16af61db-c470-4634-93b5-67d94ef6983d
# ╟─c8fd5b7e-6dce-44a5-a8a6-6c1c84dba0d2
# ╠═e9f47028-8d85-4d46-b12b-d1bff9489275
# ╟─73a19514-66f8-4508-b561-413592380098
# ╟─6244bbdd-bf0c-4d2f-bf8d-ba5b6c0daf1a
# ╟─5aa6c145-f149-49ae-8644-64d5d6ae7ac1
# ╠═506d9c37-3e08-477c-a6b8-dc735c271672
# ╟─9ee859b8-84b4-4ff9-84a9-6c20244313bb
# ╟─1ec7f614-aabd-46d2-b392-c3c3f2c5093a
# ╟─0e51622e-0f0b-41db-acee-487ea486d75f
# ╠═1b31714e-6b07-4cee-9205-a907a4fae9ef
# ╟─4d08803b-9e83-42f7-aec9-f2b4d0fca661
# ╟─e2b59fe6-0a6a-4980-b390-79699db977ca
# ╟─2eb2a0c0-fdae-4290-994b-8dc39b9f88ca
# ╠═eb907c55-6357-49b9-b022-542a24741585
# ╟─2211cd3c-1a8b-49b7-ad39-1c7849a5875e
# ╟─eb4669a5-aecb-4320-8477-9ef581a6db2f
# ╠═1f568201-92d4-493e-a9d0-7d8ca42caeb9
# ╟─8ce20d39-2a0c-490c-97f9-2e21559a6bdf
# ╟─3dd054b6-65a8-48c5-a91c-14010aa06725
# ╟─ff459b9c-41c4-4ef4-802c-fd486b5e1504
# ╠═24dc0315-e841-4eb6-8d3b-b8c6ba5c91f8
# ╟─76dcf498-0fe5-467c-bbb7-ed66d34e7f46
# ╟─a833b314-b7a3-4bf4-a092-deb9875d3a7b
# ╟─f92403c6-75bf-4886-bd96-96cd8aa94146
# ╠═f3650b49-a046-4caf-b97f-51a476764fff
# ╟─02130696-3615-42c2-a4c0-5d15a54d2f3a
# ╟─b824b498-c3ec-481a-89d7-922bcc9e273f
# ╟─4b20e231-420b-4813-ae9e-8f439773ed78
# ╟─23fa569d-6fb4-4dd8-9a62-e7f2abd85de3
# ╠═9ef995ba-cd4a-4ab0-98d8-42a68f832dd6
# ╟─83619cb0-6edb-46ff-8367-0378b2b2783e
# ╟─da26dc03-a2c2-4333-a0ee-403039d20c65
# ╟─4dbade62-2d8d-46b0-83ad-1c3982ddf5c9
# ╟─7db50062-1c29-4768-a9f5-7a24259850e3
# ╟─7e657d06-8fa1-4fc2-8427-3293701bff60
# ╟─2049e27b-bbc0-4e33-9e2c-f1f0842e97f3
# ╟─0c64314d-3c2b-43df-85e1-121e3ad2c2ba
# ╟─a6d9cb7c-28f1-48ad-8ff3-aaa872ef5da8
# ╟─611fbf7c-51de-457f-a35f-230db8992979
# ╟─d1848d82-a027-4e1d-bdff-9a60a8a1986b
# ╟─9f9ab022-e122-4880-b2a5-419def02de86
# ╟─7dc65fe5-575f-458f-a562-26f0adc380c0
# ╟─aed8567f-b982-4da6-a913-69091fe49962
# ╟─bf9fd12b-6361-4188-8bf4-1db0fc001679
# ╟─db154864-5c13-46f4-8460-3764ec871f40
# ╟─49a21c62-8ff0-4a80-93e0-b82f801f6d00
# ╟─513808b0-bb41-40b4-b805-a2e33dcaf8b5
# ╟─11df28d5-73a7-4bd0-9303-baf84f916668
# ╟─43a0bf36-ed3d-49af-9be2-a26295a977f6
# ╟─6e9f815b-fd56-4bd3-a2d0-ccc89ee81c64
# ╟─d9b1484a-c91c-4def-956c-3fb05469c6ce
# ╟─d524e5fb-0622-4503-8951-dbfc2be84088
# ╟─8ba7a253-482e-44f8-b7e3-5a3485bf343e
# ╟─f2497cb5-efd6-4045-80d2-63691055a70b
# ╟─e6d8b403-f8d0-444e-a994-e927d54b7d80
# ╟─4b4da494-897c-4a21-92eb-65f59901fec0
# ╟─37c19400-6f85-4141-ba1b-d70e7679a2be
# ╟─60fd07be-18fe-43bb-934a-8467a6753fd5
# ╟─54bbc506-0c18-4ab3-a3c7-1d811c54f73e
# ╟─5a600020-6012-4ba4-8ed2-cd08ec7bc4d5
# ╟─13f694d0-9ea8-413c-943a-b196baad5d89
# ╟─39b151a7-7d0b-4474-9c94-289136fe4c50
# ╟─07644fb0-3db1-4e46-adf6-5bac7cdb3bf3
# ╟─d7951d3a-2d7b-4303-aa3d-0b326ea51e53
# ╟─38fe75d7-ee09-4dbc-aedd-fbf8ef304bb9
# ╟─284edea4-bb09-4efb-9f5a-de9bda495914
# ╟─db1b70dc-9f0c-473d-94fc-2d87dacc1bc9
# ╟─1b039417-c067-4ec6-9047-7b1f27ca16ef
# ╠═58498690-62ba-46e9-955c-f964970cb052
# ╟─60138269-4108-44a4-9b4d-dcf397931be0
# ╟─56acf259-9998-4d39-a99b-2bf4cfaf37a6
# ╟─f461a848-0d30-4fb6-81ca-27905544d4e0
# ╟─5b5243b0-189f-48ca-a8c2-5aa63f7c74c1
# ╟─87a57724-8403-44ed-85ef-6d101409b299
# ╟─61b38b1c-e9be-43f6-a642-be6b72a27d5c
# ╠═50e0e483-82f9-4b58-baac-ea0df32b8cfa
# ╟─5f9d2291-2d6b-4a0a-a19d-c3a9ac28707e
# ╟─de186f85-91b3-4dd8-9f66-aca0ade1cb03
