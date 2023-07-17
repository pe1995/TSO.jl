### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e6969ec0-11b3-11ee-397f-dd62d15b32cb
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); Pkg.add("LaTeXStrings")
	using PythonPlot
	using TSO
	using LaTeXStrings

	plt = pyplot;
end

# ╔═╡ 60354562-1c5a-43ec-bd2d-b5b848aae452
md"## Load tables"

# ╔═╡ fd05c55a-32a8-44ab-82d9-d3da150b8228
table = abspath("tables/TSO_MARCS_v1.6")

# ╔═╡ 5dff65de-20fc-495d-b65f-1f8873c2cf8f
binned_table = abspath("DIS_MARCS_v1.6.8")

# ╔═╡ d9e4f097-6ae7-4d82-a8f5-9653be33ab7a
eos_raw = reload(
	SqEoS, joinpath(table, "combined_ross_eos.hdf5")
)

# ╔═╡ 7ef8975a-3c95-41cf-8cb1-1ad76d8d2595
opa_raw = reload(
	SqOpacity, joinpath(table, "combined_opacities.hdf5"), mmap=true
)

# ╔═╡ 96f0111b-1b07-4c4a-b92c-f2adecc809f9
size(opa_raw.κ)

# ╔═╡ fbc54d19-250a-4696-bf7a-cb3b48a680bf
eos = reload(
	SqEoS, joinpath(binned_table, "eos.hdf5")
)

# ╔═╡ 4d1f672c-b9a7-4702-9a19-33718d0a0e2e
opa = reload(
	SqOpacity, joinpath(binned_table, "binned_opacities.hdf5")
)

# ╔═╡ ad97f52a-8d8e-4991-af93-4e19c824f0a2
md"## Model
For the opatical depth of the model we use the unbinned table."

# ╔═╡ 0129fe06-20b5-4c85-bb15-b07cd731e413
# ╠═╡ show_logs = false
solar_model = upsample(
	@optical(Average3D(eos_raw, "sun_stagger.dat"), eos_raw, opa_raw), 2000
)	

# ╔═╡ 40facadb-ca4e-4b3b-b2de-754fa38b385a
md"# Radiative transfer
Compute binned and unbinned radiative transfer and compare the results."

# ╔═╡ f4e99fb7-d33e-4540-9c5a-d8192ae5a75f
weights = ω_midpoint(opa_raw)

# ╔═╡ 4d9f224f-22a3-424c-b204-c8f2af42f713
md"Applying that the opa opacities are in fact binned"

# ╔═╡ 82244213-1ceb-4e20-9ad6-cdc5ac85c259
opacities = @binned opa

# ╔═╡ 90466f84-6bb7-4fd7-978e-a86c8a4e5482
md"Applying that the opa_raw opacities are in fact unbinned"

# ╔═╡ bc7f8d42-caec-41a9-a3a2-d45e7902ec7b
opacities_raw = @binned opa_raw eos_raw

# ╔═╡ 43205993-2ee9-4677-980a-25852229eb99
solver_raw = Solver(solar_model, @axed(eos_raw), opacities=opacities_raw)

# ╔═╡ fcf83b50-f378-451d-bfe1-c4b460bfd2de
solver = Solver(solar_model, @axed(eos), opacities=opacities)

# ╔═╡ 09a5c20e-c790-4c0a-9f47-224461c3e41a
md"Solve the radiative transfer using those solvers"

# ╔═╡ e078b109-79af-44bd-80f3-ec29083144d5
@show typeof(solver)

# ╔═╡ 67f854f5-8afb-40ba-9d94-7339633e3843
q = Qr(solver) 

# ╔═╡ b28222dd-0d9c-4eb6-af43-5dff80267af7
md"Additionally pass the weigts to the raw solver, because the integration over λ is not perfomed yet in the unbinned case."

# ╔═╡ daef203a-603c-437c-82e1-bdf3f59c1f5c
q_raw = Qr(solver_raw, weights) 

# ╔═╡ 85937574-dc04-4c8d-8c4f-828f3f587c03
z, lnT, τ = solver_raw.model[:, 1], solver_raw.model[:, 2], reverse(solar_model.τ)

# ╔═╡ ae1d4866-355e-4805-b643-56937d0372e5
begin
	plt.close()
	
	ff, axf = plt.subplots(figsize=(6,6))

	mask = log10.(τ) .< 5
	
	axf.plot(
		log10.(τ[mask]), q[mask], 
		label="binned", color="k"
	)
	
	axf.plot(
		log10.(τ[mask]), q_raw[mask], 
		label="unbinned", marker=".", ls=""
	)
	
	axf.set_ylabel(L"\rm qr", fontsize="large")
	axf.set_xlabel(L"\rm \log\ \tau_{ross}", fontsize="large")
	
	axf.minorticks_on()
	axf.tick_params(top=true, right=true, direction="in", which="both")

	axf.legend(framealpha=0, loc="lower left", fontsize="large")
	
	ff.savefig(
		joinpath(binned_table, "binning_evaluation_qr.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ 1c8b521a-331d-41d5-af51-843c0bf41a27
begin
	plt.close()

	fr, axr = plt.subplots(figsize=(6,6))

	axr.plot(
		log10.(τ[mask]), (q[mask] - q_raw[mask]) ./q_raw[mask], 
		color="k", marker="x", ls="-", markersize=5
	)
	
	axr.set_ylabel(L"\rm \left(qr_{bin} - qr\right)\ /\ qr", fontsize="large")
	axr.set_xlabel(L"\rm \log\ \tau_{ross}", fontsize="large")
	
	axr.minorticks_on()
	axr.tick_params(top=true, right=true, direction="in", which="both")

	axr.set_ylim(-0.5, 0.5)

	axr.axvline(-4, color="k", ls=":", alpha=0.5)
	
	fr.savefig(
		joinpath(binned_table, "binning_evaluation_relqr.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═e6969ec0-11b3-11ee-397f-dd62d15b32cb
# ╟─60354562-1c5a-43ec-bd2d-b5b848aae452
# ╠═fd05c55a-32a8-44ab-82d9-d3da150b8228
# ╠═5dff65de-20fc-495d-b65f-1f8873c2cf8f
# ╠═d9e4f097-6ae7-4d82-a8f5-9653be33ab7a
# ╠═7ef8975a-3c95-41cf-8cb1-1ad76d8d2595
# ╠═96f0111b-1b07-4c4a-b92c-f2adecc809f9
# ╠═fbc54d19-250a-4696-bf7a-cb3b48a680bf
# ╠═4d1f672c-b9a7-4702-9a19-33718d0a0e2e
# ╟─ad97f52a-8d8e-4991-af93-4e19c824f0a2
# ╠═0129fe06-20b5-4c85-bb15-b07cd731e413
# ╟─40facadb-ca4e-4b3b-b2de-754fa38b385a
# ╠═f4e99fb7-d33e-4540-9c5a-d8192ae5a75f
# ╟─4d9f224f-22a3-424c-b204-c8f2af42f713
# ╠═82244213-1ceb-4e20-9ad6-cdc5ac85c259
# ╟─90466f84-6bb7-4fd7-978e-a86c8a4e5482
# ╠═bc7f8d42-caec-41a9-a3a2-d45e7902ec7b
# ╠═43205993-2ee9-4677-980a-25852229eb99
# ╠═fcf83b50-f378-451d-bfe1-c4b460bfd2de
# ╟─09a5c20e-c790-4c0a-9f47-224461c3e41a
# ╠═e078b109-79af-44bd-80f3-ec29083144d5
# ╠═67f854f5-8afb-40ba-9d94-7339633e3843
# ╟─b28222dd-0d9c-4eb6-af43-5dff80267af7
# ╠═daef203a-603c-437c-82e1-bdf3f59c1f5c
# ╠═85937574-dc04-4c8d-8c4f-828f3f587c03
# ╟─ae1d4866-355e-4805-b643-56937d0372e5
# ╠═1c8b521a-331d-41d5-af51-843c0bf41a27
