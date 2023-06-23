### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ cfa7c282-10e3-11ee-22a7-87824526e9a0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
	using Plots
end

# ╔═╡ 7b15f02d-a24c-43ad-9fe7-249b343f1f3f
md"# Issues with Magg et al. 2022 tables"

# ╔═╡ 7152b717-92a3-4426-839f-125459f6705d
Plots.default(fontfamily = ("Courier"), titlefont=("Courier"))

# ╔═╡ 6b6e1da8-09f3-462f-b573-df3de8c7d662
md"## The interpolated tables"

# ╔═╡ ec1b6727-97e8-4441-9a7d-c237c38c9981
magg_path = abspath("../../../tests/tables/TSO_MARCS_v1.6")

# ╔═╡ 0bf842ec-e5be-4fc1-8cd7-13e47959c418
norm_path = abspath("../../../tests/tables/TSO_MARCS_v1.4")

# ╔═╡ 7bf1afdd-3233-4e31-860d-896bfe81850f
begin
	eos_magg = reload(SqEoS, joinpath(magg_path, "combined_eos.hdf5"))
	opa_magg = reload(
		SqOpacity, joinpath(magg_path, "combined_opacities.hdf5"), mmap=true
	)
end

# ╔═╡ 62af4c18-6c2a-4780-b172-c64d9513436c
begin
	eos_norm = reload(SqEoS, joinpath(norm_path, "combined_eos.hdf5"))
	opa_norm = reload(
		SqOpacity, joinpath(norm_path, "combined_opacities.hdf5"), mmap=true
	)
end

# ╔═╡ 68d82824-0c35-40cd-bc9d-1b9b632309bb
solar_model = Average3D("../../../tests/stagger_av.dat")

# ╔═╡ 0a3fa85e-6c1d-4896-a24c-03044a7c0c66
ilam = [10:1000:100000...]

# ╔═╡ 178167d2-d888-4d40-80b8-5b99c8041810
# ╠═╡ show_logs = false
k_magg = [
	lookup(eos_magg, opa_magg, :κ, solar_model.lnρ, solar_model.lnT, i) 
	for i in ilam
]

# ╔═╡ 0224f1ce-0e9a-4ae1-90a8-0f6d21ec6a5c
k_norm = [
	lookup(eos_norm, opa_norm, :κ, solar_model.lnρ, solar_model.lnT, i-1) 
	for i in ilam
]

# ╔═╡ c2e51c3b-b5a7-4e90-8252-0e965d90b633
begin
	plot(framestyle=:box, grid=false, 
		legendforegroundcolor=nothing,
		legendbackgroundcolor=nothing)

	#plot!(title="λ = $(opa_magg.λ[ilam]) Å")

	for i in eachindex(ilam)
		if i ==1
			plot!(solar_model.lnT, log10.(k_norm[i]), 
				label="original table", color=:red)
			
			plot!(solar_model.lnT, log10.(k_magg[i]), 
				label="Magg et al. 2022", color=:black)
		else
			plot!(solar_model.lnT, log10.(k_norm[i]), 
				label=nothing, color=:red)
			
			plot!(solar_model.lnT, log10.(k_magg[i]), 
				label=nothing, color=:black)
		end
	end

	plot!(xlabel="ln T", ylabel="log κ")
	plot!(ylim=(-5, 8))
end

# ╔═╡ de7a88b3-5bb6-47b3-a17c-a4ca30920848
md"## Formation Opacities"

# ╔═╡ a32ea314-ab43-4306-b536-7cb4cb8bdf57
fopa_magg = reload(
	SqOpacity, joinpath(magg_path, "combined_formation_opacities.hdf5"),
	mmap=true
)

# ╔═╡ ce1b2e50-8c82-4a69-af84-dd52163a5f96
fopa_norm = reload(
	SqOpacity, joinpath(norm_path, "combined_formation_opacities.hdf5"),
	mmap=true
)

# ╔═╡ 590629aa-6dd6-4fd2-800c-a7d86761f0e7
begin
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing,
		legendbackgroundcolor=nothing)

	scatter!(log10.(fopa_norm.λ), -log10.(fopa_norm.κ_ross), 
		markersize=1, color=:red, label="original table",
		markerstrokecolor=:red)

	scatter!(log10.(fopa_magg.λ), -log10.(fopa_magg.κ_ross), 
		markersize=1, color=:black, label="Magg et al. 2022",
		markerstrokecolor=:black)

	plot!(ylabel="- τ-ross (τλ = 1)", xlabel="wavelength")
end

# ╔═╡ Cell order:
# ╟─7b15f02d-a24c-43ad-9fe7-249b343f1f3f
# ╠═cfa7c282-10e3-11ee-22a7-87824526e9a0
# ╠═7152b717-92a3-4426-839f-125459f6705d
# ╠═6b6e1da8-09f3-462f-b573-df3de8c7d662
# ╠═ec1b6727-97e8-4441-9a7d-c237c38c9981
# ╠═0bf842ec-e5be-4fc1-8cd7-13e47959c418
# ╠═7bf1afdd-3233-4e31-860d-896bfe81850f
# ╠═62af4c18-6c2a-4780-b172-c64d9513436c
# ╠═68d82824-0c35-40cd-bc9d-1b9b632309bb
# ╠═0a3fa85e-6c1d-4896-a24c-03044a7c0c66
# ╠═178167d2-d888-4d40-80b8-5b99c8041810
# ╠═0224f1ce-0e9a-4ae1-90a8-0f6d21ec6a5c
# ╟─c2e51c3b-b5a7-4e90-8252-0e965d90b633
# ╟─de7a88b3-5bb6-47b3-a17c-a4ca30920848
# ╠═a32ea314-ab43-4306-b536-7cb4cb8bdf57
# ╠═ce1b2e50-8c82-4a69-af84-dd52163a5f96
# ╟─590629aa-6dd6-4fd2-800c-a7d86761f0e7
