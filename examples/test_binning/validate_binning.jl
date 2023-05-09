### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ c239d5b6-ea86-11ed-3da9-436677edfc25
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
	using Plots
	using DelimitedFiles
end

# ╔═╡ ac3320db-ff4d-42d7-93de-721d12026088
md"
# Validate Binned Opacities
"

# ╔═╡ 2c39a006-5385-48cc-8ab6-cbdf10adfef0
begin
	Plots.default(fontfamily = ("Courier"), titlefont=("Courier"))

	function copy_ticks(sp::Plots.Subplot; minorticks=10)
	    ptx = twinx(sp)
	    plot!(ptx,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
		
	    pty = twiny(sp)
	    plot!(pty,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
	end
	
	copy_ticks(plt::Plots.Plot = current(); minorticks=10) = copy_ticks(plt[1], minorticks=minorticks)
	
	
	function basic_plot!(plot::Plots.Plot=current(); minorticks=10, 
						tickfontsize=12, 
						legendfontsize=12, 
						guidefontsize=12, 
						size=(1000,500), 
						lm=5, bm=5, tm=5, rm=5)
		
	    copy_ticks(plot, minorticks=15)
	    plot!(plot, framestyle=:box, 
				minorticks=minorticks, 
				tickfontsize=tickfontsize, 
				legendfontsize=legendfontsize, 
	            grid=false, 
				foreground_color_legend=nothing,
				background_color_legend=nothing,
				size=size, 
				leftmargin=lm*Plots.mm, rightmargin=rm*Plots.mm, 
				bottommargin=bm*Plots.mm, topmargin=tm*Plots.mm, 
				guidefontsize=guidefontsize)
	end;
end

# ╔═╡ 1b529b15-28e6-48b8-afa0-01f7b6fa0516
md"
## Loading Tables
We validate the binning by comparing the 1D heating rate with the unbinned solution.
First, we load the unbinned opacity table + EoS everything else will be compared to.
"

# ╔═╡ a2ad12b7-ce18-45c9-83d1-280d3867e1d6
eos_raw = reload(SqEoS, 
		joinpath("tables/TSO_MARCS_v1.4", "combined_ross_eos.hdf5"))

# ╔═╡ 7ba372f7-a686-40e7-9c52-bac48787c6f1
opa_raw = reload(SqOpacity, 
		joinpath("tables/TSO_MARCS_v1.4", "combined_opacities.hdf5"), mmap=true)

# ╔═╡ 043a9791-05e1-44ac-b325-54e281650f29
md"We have to option to compare different binnings with each other."

# ╔═╡ 08898890-d277-4d05-b10e-84798ba7fc87
folders = ["DIS_MARCS_v1.4.1"]

# ╔═╡ 0419df91-6a85-4f6e-975a-e3aac8760057
labels = ["5 bins (no λ splits)"]

# ╔═╡ 968d8c76-d0d0-4b90-a5c7-0c2a8990f2e0
eos = [reload(SqEoS, joinpath(f, "eos.hdf5")) for f in folders]

# ╔═╡ 5efed2cc-a1c4-45aa-b0d5-958b74b07b35
opa = [reload(SqOpacity, joinpath(f, "binned_opacities.hdf5")) for f in folders]

# ╔═╡ 939212b2-fcb8-4cc9-a3ce-93fc28cb3e85
md"## Solving radiative transfer"

# ╔═╡ a164a69a-e3c6-4153-a05e-3df20f43c97c
weights = ω_midpoint(opa_raw)

# ╔═╡ aea42941-06db-440a-92be-99993f1e1cfc
opacities = [@binned(opa_i) for opa_i in opa]

# ╔═╡ fd3b8e56-49e2-4ba4-999f-bc6363cb8b8d
opacities_raw = @binned opa_raw eos_raw

# ╔═╡ 36eee85a-9684-4991-8cac-01da3fb1ec8f
md"The model atmosphere we want to use for comparison. We also compute the optical depth based on the unbinned table."

# ╔═╡ 3466f47d-2c42-4281-bdeb-ec6e651e86d8
solar_model = @optical Average3D(eos_raw, "stagger_av.dat") eos_raw opa_raw;

# ╔═╡ 9057ed2d-b5f2-4fa7-bd86-539a83044f3f
md"Next we construct the solvers for these tables."

# ╔═╡ c40aaba1-0eae-4aee-9790-576f4bb4866c
solvers = [Solver(solar_model, @axed(eos[i]), opacities=opacities[i]) 
								for i in eachindex(eos)]

# ╔═╡ 3500aa35-7fb7-4d6e-9ef4-afa1d6709578
solver_raw = Solver(solar_model, @axed(eos_raw), opacities=opacities_raw)

# ╔═╡ 00c91510-716b-400c-872f-d1c391eb3cb4
md"Now solve"

# ╔═╡ 2a551262-d3a3-47fc-a90a-5503a52c8fb5
q = [Qr(solver) for solver in solvers]

# ╔═╡ 646b9ae8-348c-4550-80a6-6f7eaa3d43c6
md"In the unbinned case, we have to additionally provide λ integration weights, because the λ integration in the binned case happens during the binning already."

# ╔═╡ e9900f75-11b1-4c50-baa8-d5d2d1e9f2c6
q_raw = Qr(solver_raw, weights)

# ╔═╡ 52e2e5fc-db22-4414-9ab7-28995ee2d26c
md"# Compare Heating rates"

# ╔═╡ 524090eb-f1aa-44d6-b8f9-2699a14bd7a8
z, lnT, τ = solver_raw.model[:, 1], solver_raw.model[:, 2], reverse(solar_model.τ)

# ╔═╡ 54417cb3-a0fc-45cf-bb2a-3d413b0018d1
begin
	mask = log10.(τ) .< 5
	plot(log10.(τ)[mask], q_raw[mask], label="unbinned", lw=3, 
		ylabel="heating rate", xlabel="log τ ross")

	for i in eachindex(q)
		plot!(log10.(τ), q[i], label=labels[i], lw=3, ls=:dash)
	end

	plot!(xlim=(-8, 5.5))
	basic_plot!(size=(900,600))

	plot!(legend=:bottomleft)
end

# ╔═╡ abec9170-f0f7-42ff-a16e-126733731db1
begin
	plot(ylabel="heating rate (rel. difference)", xlabel="log τ ross")
	
	for i in eachindex(q)
		plot!(log10.(τ[mask]), (q[i][mask] .- q_raw[mask]) ./ q_raw[mask] , 
				label=labels[i], lw=3)
	end

	plot!(xlim=(-8, 5.5))
	basic_plot!(size=(900,600))

	plot!(legend=:topright)
end

# ╔═╡ Cell order:
# ╟─ac3320db-ff4d-42d7-93de-721d12026088
# ╠═c239d5b6-ea86-11ed-3da9-436677edfc25
# ╟─2c39a006-5385-48cc-8ab6-cbdf10adfef0
# ╟─1b529b15-28e6-48b8-afa0-01f7b6fa0516
# ╠═a2ad12b7-ce18-45c9-83d1-280d3867e1d6
# ╠═7ba372f7-a686-40e7-9c52-bac48787c6f1
# ╟─043a9791-05e1-44ac-b325-54e281650f29
# ╠═08898890-d277-4d05-b10e-84798ba7fc87
# ╠═0419df91-6a85-4f6e-975a-e3aac8760057
# ╠═968d8c76-d0d0-4b90-a5c7-0c2a8990f2e0
# ╠═5efed2cc-a1c4-45aa-b0d5-958b74b07b35
# ╟─939212b2-fcb8-4cc9-a3ce-93fc28cb3e85
# ╠═a164a69a-e3c6-4153-a05e-3df20f43c97c
# ╠═aea42941-06db-440a-92be-99993f1e1cfc
# ╠═fd3b8e56-49e2-4ba4-999f-bc6363cb8b8d
# ╟─36eee85a-9684-4991-8cac-01da3fb1ec8f
# ╠═3466f47d-2c42-4281-bdeb-ec6e651e86d8
# ╟─9057ed2d-b5f2-4fa7-bd86-539a83044f3f
# ╠═c40aaba1-0eae-4aee-9790-576f4bb4866c
# ╠═3500aa35-7fb7-4d6e-9ef4-afa1d6709578
# ╟─00c91510-716b-400c-872f-d1c391eb3cb4
# ╠═2a551262-d3a3-47fc-a90a-5503a52c8fb5
# ╟─646b9ae8-348c-4550-80a6-6f7eaa3d43c6
# ╠═e9900f75-11b1-4c50-baa8-d5d2d1e9f2c6
# ╟─52e2e5fc-db22-4414-9ab7-28995ee2d26c
# ╠═524090eb-f1aa-44d6-b8f9-2699a14bd7a8
# ╟─54417cb3-a0fc-45cf-bb2a-3d413b0018d1
# ╟─abec9170-f0f7-42ff-a16e-126733731db1
