### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 2c8f86fe-68f5-11ee-1788-9f6f32429204
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using PythonPlot
	using TSO
	using LaTeXStrings

	plt = pyplot
end;

# ╔═╡ 607818f0-d5a7-4644-8045-36a2c8f396cf
matplotlib.style.use(joinpath(dirname(pathof(TSO)), "Bergemann2023.mplstyle"))

# ╔═╡ c99d2f52-9c05-4492-81b3-64b61ab52068
md"## Load tables"

# ╔═╡ 0ab293b6-8641-4148-aaae-6fe76d280eec
table = abspath("tables/TSO_MARCS_v1.6")

# ╔═╡ eb112e63-7440-4b94-97cf-f2f87d0cc4ac
name = "t50g40m00"
#name = "cnt"

# ╔═╡ b86e8a40-6d20-4511-b0b2-75bc0837fae9
binned_table = abspath("../../../../model_grid/examples/initial_models/DIS_MARCS_$(name)_v0.1")
#binned_table = abspath("DIS_MARCS_cnt_v1.7.6") 

# ╔═╡ a3f3b014-8358-48c0-88ff-185a12525685
eos_raw = TSO.reload(
	TSO.SqEoS, joinpath(table, "combined_ross_eos.hdf5")
)

# ╔═╡ fd60893c-fdcf-4ffd-8ec2-89acef5c82ee
opa_raw = TSO.reload(
	TSO.SqOpacity, joinpath(table, "combined_opacities.hdf5"), mmap=true
)

# ╔═╡ da4b3c11-a168-4823-b8ba-7ab639251b04
size(opa_raw.κ)

# ╔═╡ 33720dc6-cf40-4be5-aafc-41b51280d4ed
eos = TSO.reload(
	TSO.SqEoS, joinpath(binned_table, "eos.hdf5")
)

# ╔═╡ 546500e8-4d8d-46cc-9d2b-8e8590bf8310
opa = TSO.reload(
	TSO.SqOpacity, joinpath(binned_table, "binned_opacities.hdf5")
)

# ╔═╡ 5fcb0bb7-88ac-4e03-8357-63d7d4ec9cd9
size(opa.κ)

# ╔═╡ 6de74388-7d45-48e4-8982-167954b6634e
md"## Model
For the opatical depth of the model we use the unbinned table."

# ╔═╡ 88bc870e-a6b4-4630-9ef5-d14aaf316127
function cut(m::TSO.AbstractModel; kwargs...)
	masks = trues(length(m.z))
	for (para, lims) in kwargs
		v = getfield(m, para)
		mask = first(lims) .< v .< last(lims)
		masks .= masks .& mask
	end

	dat = []
	for f in fieldnames(typeof(m))
		v = getfield(m, f)
		if typeof(v) <: AbstractArray
			append!(dat, [v[masks]])
		else
			append!(dat, [v])
		end
	end

	typeof(m)(dat...)
end

# ╔═╡ eadda4ab-e571-4e20-aef2-f8b82f376626
solar_model_up = TSO.upsample(
	TSO.@optical(TSO.Average3D(eos_raw, "sun_stagger.dat"), eos_raw, opa_raw), 2000
)

# ╔═╡ 0ab8865b-15a5-4c90-a69d-57c287d69d9e
solar_model_orig = TSO.@optical(
	TSO.Average3D(eos_raw, "sun_stagger.dat"), eos_raw, opa_raw
)

# ╔═╡ 6cbe21bb-293c-4f2e-b6a6-6c02dc797a21
solar_model = solar_model_orig

# ╔═╡ f736ebf6-3e3a-49b5-b14b-5624c272b0c0
solar_model.z

# ╔═╡ 1e9365c2-26b1-4075-a560-5b5a09d020ca
md"# Radiative transfer
Compute binned and unbinned radiative transfer and compare the results."

# ╔═╡ a9c160d6-673f-4cfa-8358-4c01be07d1e6
weights = TSO.ω_midpoint(opa_raw)

# ╔═╡ 243a7aee-cd48-4e41-ad46-9e255019e020
md"Applying that the opa opacities are in fact binned"

# ╔═╡ a534574a-643d-47c0-899a-1be01aa2e5cb
opacities = TSO.@binned opa

# ╔═╡ b04293c5-0ed9-486c-831b-4ed4bd62f31c
md"Applying that the opa_raw opacities are in fact unbinned"

# ╔═╡ bff0d9f0-fb5e-412d-8e58-35636a0cf888
opacities_raw = TSO.@binned opa_raw eos_raw

# ╔═╡ 84d31a8e-d586-4422-b636-74ec0111b608
md"The solvers"

# ╔═╡ 36ea0356-223a-4b4c-82ca-a755800a8e03
solver = TSO.Heating(solar_model, eos=TSO.@axed(eos), opacities=opacities)

# ╔═╡ 4369fb89-fca6-46ad-bbb4-5eae583e6f43
solver_raw = TSO.Heating(solar_model, eos=TSO.@axed(eos_raw), opacities=opacities_raw)

# ╔═╡ e71fac10-3069-457c-be96-d64d3ec4ba10
md"# Heating"

# ╔═╡ 33fb636a-e216-425b-915d-4a87c0ed62ce
Q = TSO.solve(solver)

# ╔═╡ 7f76414d-4db8-4d2b-9f59-0aacf6635e8d
Q_raw = TSO.solve(solver_raw, weights=weights)

# ╔═╡ 9790bf4a-f13c-4c89-9f3a-139080ed4dbd
z, lnT, τ = Q.model.z, Q.model.lnT, Q.model.τ

# ╔═╡ fb167940-25a4-4454-8cea-68e8c32fc19c
md"# Comparison"

# ╔═╡ 4b33966c-d8ed-48dc-b520-93c4c4b199c0
begin
	plt.close()
	
	ff, axf = plt.subplots(figsize=(6,6))

	mask = -5 .< log10.(τ) .< 5

	# overplot the interpolation lines
	#=for line in solar_model_orig.τ
		axf.axvline(log10.(line), alpha=0.3)
	end=#

	axf.plot(
		log10.(τ[mask]), 
		#z[mask],
		TSO.heating(Q_raw)[mask] ,#./ exp.(TSO.model(Q_raw).lnρ[mask]), 
		label="unbinned", marker=".", ls="", markersize=8
	)

	axf.plot(
		log10.(τ[mask]), 
		#z[mask],
		TSO.heating(Q)[mask] ,#./ exp.(TSO.model(Q).lnρ[mask]), 
		label="binned", color="k"
	)
	
	axf.set_ylabel(L"\rm qr", fontsize="large")
	axf.set_xlabel(L"\rm \log\ \tau_{ross}", fontsize="large")
	
	axf.minorticks_on()
	axf.tick_params(top=true, right=true, direction="in", which="both")

	axf.set_xlim(-4, 4)
	#axf.set_yscale("log")
	
	axf.legend(framealpha=0, loc="lower left", fontsize="large")
	
	ff.savefig(
		joinpath(binned_table, "binning_evaluation_qr.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ 11be2ede-35f5-445f-bdcf-24b281cf1871
begin
	plt.close()

	fr, axr = plt.subplots(figsize=(5,6))

	xr = TSO.heating(Q)[mask] #./ exp.(TSO.model(Q).lnρ[mask])
	yr = TSO.heating(Q_raw)[mask] #./ exp.(TSO.model(Q_raw).lnρ[mask])
	m = sortperm(xr)

	qr_max = maximum(abs.(yr))
	xr = xr ./ qr_max
	yr = yr ./ qr_max
	
	axr.plot(
		range(minimum(xr), maximum(xr), length=500) |> collect,
		range(minimum(xr), maximum(xr), length=500) |> collect,
		color="k", ls=":"
	)
	
	axr.plot(
		xr[m], yr[m],
		color="k", marker="s", ls="", markersize=6
	)
	
	axr.set_ylabel(
		L"\rm q_{rad}^{\lambda}\ /\ \left| q_{rad}^{\lambda,max} \right|", 
		fontsize="large"
	)
	axr.set_xlabel(
		L"\rm q_{rad}^{bin}\ /\ \left| q_{rad}^{\lambda,max} \right|",
		fontsize="large"
	)
	
	axr.minorticks_on()
	axr.tick_params(top=true, right=true, direction="in", which="both")
	
	fr.savefig(
		joinpath(binned_table, "binning_evaluation_qqr.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ 65ec74da-3809-4194-90e8-ca0bcc20adc6
begin
	plt.close()

	fd, axd = plt.subplots(figsize=(6,6))

	xd = TSO.heating(Q)[mask] #./ exp.(TSO.model(Q).lnρ[mask])
	yd = TSO.heating(Q_raw)[mask] #./ exp.(TSO.model(Q_raw).lnρ[mask])

	axd.plot(
		log10.(τ[mask]),
		(xd .- yd) ./ maximum(abs.(yd)) .* 100,
		color="k", ls="", marker="s", markersize=7, 
		markeredgecolor="k", markerfacecolor="None"
	)
	
	axd.set_ylabel(
		L"\rm \left( q_{rad}^{bin} - q_{rad}^{\lambda} \right)\ /\ \left| q_{rad}^{\lambda, max} \right| \ [\%]", 
		fontsize="large"
	)
	axd.set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")
	
	axd.minorticks_on()
	axd.tick_params(top=true, right=true, direction="in", which="both")

	axd.set_ylim(-4, 4)
	
	fd.savefig(
		joinpath(binned_table, "binning_evaluation_relqr.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ ba1497d9-3c4e-4bae-a7d5-46f8eaf956f8
begin
	plt.close()

	fs, axs = plt.subplots(figsize=(6,6))

	xs = TSO.heating(Q)[mask] #./ exp.(TSO.model(Q).lnρ[mask])
	ys = TSO.heating(Q_raw)[mask] #./ exp.(TSO.model(Q_raw).lnρ[mask])
	
	ss = (xs .- ys) ./ ys
	ms = -20 .< ss .< 20

	means = round(TSO.mean(ss[ms]), sigdigits=3)
	stds = round(TSO.std(ss[ms]), sigdigits=2)
	maxd = round(maximum(abs.(xs .- ys)) ./ maximum(abs.(ys)), sigdigits=3)
	
	axs.hist(
		ss[ms],
		label="mean: $(means)\nstd: $(stds)\nmax dq: $(maxd)",
		color="k", bins=100
	)
	
	axs.set_xlabel(L"\rm qr_{bin}\ /\ q - 1", fontsize="large")
	
	axs.minorticks_on()
	axs.tick_params(top=true, right=true, direction="in", which="both")

	axs.set_xlim(-.5, .5)
	axs.legend(framealpha=0)
	
	fs.savefig(
		joinpath(binned_table, "binning_evaluation_stat.png"), 
		dpi=800, bbox_inches="tight"
	)
	
	gcf()
end

# ╔═╡ fdac8f67-3712-4a48-9041-c07f92fabc19
md"# Assignment"

# ╔═╡ 4dcf87c7-ef23-4f13-b158-dc8436937fe7
begin
	# the formation opacities are stored in the unbinned table
	fopaD = TSO.reload(
		TSO.SqOpacity, 
		joinpath(
			table, "combined_formation_opacities_$(name)_magg22.hdf5"
		),
		mmap=true
	)
	

	# the bin assignment should be saved for this opacity table
	fid = TSO.HDF5.h5open(
		joinpath(binned_table, "bin_assignment.hdf5"), 
		"r"
	)
	
	binassignment = TSO.HDF5.read(fid["bins"])
	λ = TSO.HDF5.read(fid["lambda"])
	close(fid)


	
	# Figure
	cD = plt.cm.gnuplot2
	normD = matplotlib.colors.BoundaryNorm(
		range(0.5, length(opa.λ)+0.5, step=1) |> collect, 
		cD.N - 50
	)
	
	plt.close()
	fD, axD = plt.subplots(1, 1, figsize=(5,6))
	#visual.basic_plot!(axD)
	
	im = axD.scatter(
		log10.(λ), -log10.(fopaD.κ_ross), 
		c=binassignment,
		s=1,
		cmap=cD,
		norm=normD,
	)
	cbarD = fD.colorbar(
		im, 
		ax=axD, 
		ticks=range(1, length(opa.λ), step=1) |> collect
	)


	axD.set_ylabel(
		L"\rm - \log \tau_{ross}\ (\log \tau_{\lambda}=1)"
	)
	axD.set_xlabel(
		L"\rm \log \lambda\ [\AA]"
	)
	cbarD.set_label(
		L"\rm opacity\ bins"
	)
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═2c8f86fe-68f5-11ee-1788-9f6f32429204
# ╠═607818f0-d5a7-4644-8045-36a2c8f396cf
# ╟─c99d2f52-9c05-4492-81b3-64b61ab52068
# ╠═0ab293b6-8641-4148-aaae-6fe76d280eec
# ╠═eb112e63-7440-4b94-97cf-f2f87d0cc4ac
# ╠═b86e8a40-6d20-4511-b0b2-75bc0837fae9
# ╠═a3f3b014-8358-48c0-88ff-185a12525685
# ╠═fd60893c-fdcf-4ffd-8ec2-89acef5c82ee
# ╠═da4b3c11-a168-4823-b8ba-7ab639251b04
# ╠═33720dc6-cf40-4be5-aafc-41b51280d4ed
# ╠═546500e8-4d8d-46cc-9d2b-8e8590bf8310
# ╠═5fcb0bb7-88ac-4e03-8357-63d7d4ec9cd9
# ╟─6de74388-7d45-48e4-8982-167954b6634e
# ╟─88bc870e-a6b4-4630-9ef5-d14aaf316127
# ╠═eadda4ab-e571-4e20-aef2-f8b82f376626
# ╠═0ab8865b-15a5-4c90-a69d-57c287d69d9e
# ╠═6cbe21bb-293c-4f2e-b6a6-6c02dc797a21
# ╠═f736ebf6-3e3a-49b5-b14b-5624c272b0c0
# ╟─1e9365c2-26b1-4075-a560-5b5a09d020ca
# ╠═a9c160d6-673f-4cfa-8358-4c01be07d1e6
# ╟─243a7aee-cd48-4e41-ad46-9e255019e020
# ╠═a534574a-643d-47c0-899a-1be01aa2e5cb
# ╟─b04293c5-0ed9-486c-831b-4ed4bd62f31c
# ╠═bff0d9f0-fb5e-412d-8e58-35636a0cf888
# ╟─84d31a8e-d586-4422-b636-74ec0111b608
# ╠═36ea0356-223a-4b4c-82ca-a755800a8e03
# ╠═4369fb89-fca6-46ad-bbb4-5eae583e6f43
# ╟─e71fac10-3069-457c-be96-d64d3ec4ba10
# ╠═33fb636a-e216-425b-915d-4a87c0ed62ce
# ╠═7f76414d-4db8-4d2b-9f59-0aacf6635e8d
# ╠═9790bf4a-f13c-4c89-9f3a-139080ed4dbd
# ╟─fb167940-25a4-4454-8cea-68e8c32fc19c
# ╟─4b33966c-d8ed-48dc-b520-93c4c4b199c0
# ╟─11be2ede-35f5-445f-bdcf-24b281cf1871
# ╟─65ec74da-3809-4194-90e8-ca0bcc20adc6
# ╟─ba1497d9-3c4e-4bae-a7d5-46f8eaf956f8
# ╟─fdac8f67-3712-4a48-9041-c07f92fabc19
# ╟─4dcf87c7-ef23-4f13-b158-dc8436937fe7
