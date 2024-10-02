### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ c65e0994-8085-11ef-3cbc-fb515d8af323
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end;

# ╔═╡ cd40eac0-d41d-4c45-b7ca-d7f3cc168829
md"# MANCHA Opacities + M3DIS EoS"

# ╔═╡ 3325eff5-96ed-4823-b87e-185e75b3db97


# ╔═╡ 830393f6-65c8-49a8-bafc-447f142d04f3
path = "kappa_modelingG2V_tau_optimal4bins_60x35_synspec_trhogrid_isoROSSLENAD.h5"

# ╔═╡ be615e70-389b-47a2-9a33-6b82dda1e650
fid = TSO.h5open(path, "r")

# ╔═╡ cdcd8c91-40d2-40f6-8865-306f45e6008f
data = Dict(
	k => TSO.read(fid, k) for k in keys(fid)
)

# ╔═╡ c3d7d282-9e6e-4c8f-885e-71bdebe210de


# ╔═╡ c9ea9073-e973-42ba-b412-8a9938b11a25
begin
	@show size(data["tab_rho"])
	@show size(data["tab_t"])
	@show size(data["kap_mean"])
	@show size(data["b_band"])
	@show size(data["level"])
	@show size(data["wl_lev"])
end;

# ╔═╡ 9ae7c1f5-a2bd-4518-b686-d53b81ad9fcb


# ╔═╡ d45f9356-94f7-4320-ab75-24ef0b18fc80
md"For `TSO` we need the tables in the format (nT, n$\rm\rho$, n$\rm\lambda$), so we permute the dimensions and are ready to assign an EoS. The source function needs to be provided on the same grid, so we puff up the array to a 2D array with the same values for all densities.
If in the future an EoS will be provided, the interpolation step can be skipped."

# ╔═╡ 7469fefc-b084-4dc1-9c1a-78f747bd6eb0
begin
	Tp = Float32
	T = Base.convert.(Tp, data["tab_t"][:])
	ρ = Base.convert.(Tp, data["tab_rho"][:])
	κ = Base.convert.(Tp, permutedims(data["kap_mean"], (2, 3, 1)))
	S = Base.convert.(Tp, similar(κ))
	λ = Base.convert.(Tp, range(1., size(S, 3), step=1.) |> collect)
	κ_ross = Base.convert.(Tp,ones(length(T), length(ρ)))
	for b in axes(S, 3)
		TSO.puff_up!(view(S, :, :, b), data["b_band"][b, :])
	end
end

# ╔═╡ c4d133d7-5a9e-44b9-abcb-c544e0979135
opa = SqOpacity(κ, κ_ross, S, λ, false) 

# ╔═╡ f54bbdb9-d9e2-47b0-aafe-5302b8510290


# ╔═╡ 1ca76166-5719-470c-8cd4-c44438310c86
md"A fake EoS is needed to allow for interpolation to legit EoS."

# ╔═╡ 37e70e5a-eadd-4e3f-b6dd-e3ac05595f7e
fake_eos = SqEoS(
	log.(exp10.(ρ)), 
	log.(exp10.(T)), 
	zeros(Tp, length(T), length(ρ)), 
	zeros(Tp, length(T), length(ρ)), 
	zeros(Tp, length(T), length(ρ)),
	zeros(Tp, length(T), length(ρ))
)

# ╔═╡ 2925f271-54d7-4b83-a653-b157ac61f171


# ╔═╡ 04cc1bc5-2586-4102-914e-c1918b8868a7
md"Multiply by the table density to make sure opacity is in the right units of inverse centimeters."

# ╔═╡ 1c37cc21-b4b8-4ffb-92d9-696de4cfa462
opa_binned = TSO.@binned opa fake_eos

# ╔═╡ 7df487d7-ca05-4fa8-a7ea-62298e25d400


# ╔═╡ 1f8b2a1f-4c51-46cf-9d0e-86bf1a6dba0b
md"The opacity table above needs to be interpolated to a different EoS in case there is non provided."

# ╔═╡ 8a2f3599-c7be-4159-b56d-315b8f9b1996
eos_new_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/opacity_tables/TSO_M3D_magg_m0_a0_vmic1_v5.1/combined_eos_magg_m0_a0_vmic1.hdf5"

# ╔═╡ ff6559ef-9b86-40df-9ad6-4b921d934a6f
eos_new = reload(SqEoS, eos_new_path)

# ╔═╡ 4012a87c-03d6-45fe-b08e-41e0c3fb28a0


# ╔═╡ 27ce24a3-dca7-429e-849c-9f8accda58e4
opa_interpolated = complement(
	fake_eos, eos_new, opa_binned.opacities
)

# ╔═╡ 1ad3c07e-729c-480b-89b4-54e1f0f38f55
TSO.transfer_rosseland!(eos_new, opa_interpolated);

# ╔═╡ 42052bfb-8321-45c4-af8e-3d053e2ba623


# ╔═╡ 8d07f251-9362-495d-b502-c7b98a318140
let
	f, ax = plt.subplots(2, 1, sharex=true, sharey=true)
	plt.subplots_adjust(hspace=0.2)

	bin = 1 
	
	ax[1].set_title("M3DIS EoS")
	tt, rr = meshgrid(@axed(eos_new))
	im = ax[1].pcolormesh(
		log10.(exp.(tt)), log10.(exp.(rr)), 
		log10.(opa_interpolated.κ[:, :, bin]), 
		rasterized=true
	)
	c1 = f.colorbar(im, ax=ax[1])
	c1.set_label("opacity (bin $(bin))")
	

	ax[0].set_title("MANCHA EoS")
	tt, rr = meshgrid(@axed(fake_eos))
	im = ax[0].pcolormesh(
		log10.(exp.(tt)), log10.(exp.(rr)), 
		log10.(opa_binned.opacities.κ[:, :, bin]), 
		rasterized=true
	)
	c0 = f.colorbar(im, ax=ax[0])
	c0.set_label("opacity (bin $(bin))")

	ax[1].set_xlabel("log temperature")
	ax[1].set_ylabel("log density")
	ax[0].set_ylabel("log density")
	
	f
end

# ╔═╡ 883fcfc6-bc40-4011-ada5-031259f42320


# ╔═╡ 37630f68-acd6-4517-a8c8-58f825b00c93
eos_dir_new = "MANCHA_M3D_magg_m0_a0_vmic1_v5.1/"

# ╔═╡ d2da6ed8-d19c-4b12-ab33-384d0b524cb3
!isdir(eos_dir_new) && mkdir(eos_dir_new)

# ╔═╡ fb0dbfdf-d2ce-4b20-bdeb-fa7cb6599378
TSO.save(eos_new, joinpath(eos_dir_new, "eos.hdf5"))

# ╔═╡ 3b817455-d43d-4f61-ac3d-bdb80e992e34
TSO.save(opa_interpolated, joinpath(eos_dir_new, "binned_opacities.hdf5"))

# ╔═╡ d6ba4c07-6c86-48de-979e-e212f22da557


# ╔═╡ 445da383-5dfe-427e-89f9-28ab480ffa47
md"These tables can now be used in the usual routines, as if they were from M3D or TS. For dispatch, we need the tables on the energy grid. This is done automatically when you use this EoS (+opacities) to create a new initial model. One can however also directly do the conversion now. If one wants to make sure, that the internal energy is limited to the actual range of the model, one needs to load an average model to do so, in this case the sun."

# ╔═╡ e566b981-c20d-4405-bbd3-7f59ecdcc915
av_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/MUST.jl/initial_grids/Stagger/av_models/t5777g44m0005_00070_av.dat"

# ╔═╡ b3551f27-a745-4816-a173-6fcc54a2f001
eos_dir_new_E = "MANCHA_M3D_E_magg_m0_a0_vmic1_v5.1/"

# ╔═╡ 46f48e0e-e324-439c-9b71-762989b2b778
TSO.convert_fromT_toE(eos_dir_new, eos_dir_new_E, av_path; lnEi_stretch=1.0)

# ╔═╡ 87913ced-c8d0-42c1-916a-9250cdde8f73


# ╔═╡ 5eff0a2d-2e20-4699-a29f-d1be0cc4c428
eosE = reload(SqEoS, joinpath(eos_dir_new_E, "eos.hdf5"))

# ╔═╡ 5846c57c-171b-4ef5-abcf-666596605e6d
opaE = reload(SqOpacity, joinpath(eos_dir_new_E, "binned_opacities.hdf5"))

# ╔═╡ 2d48d193-b822-46e2-a605-397b27634b3e
let
	
	f, ax = plt.subplots(1, 1)

	bin = 1 
	
	ee, rr = meshgrid(@axed(eosE))
	im = ax.pcolormesh(
		log10.(exp.(ee)), log10.(exp.(rr)), 
		log10.(opaE.κ[:, :, bin]), 
		rasterized=true
	)
	c1 = f.colorbar(im, ax=ax)
	c1.set_label("opacity (bin $(bin))")

	ax.set_xlabel("log internal energy")
	ax.set_ylabel("log density")
	
	f

end

# ╔═╡ Cell order:
# ╟─cd40eac0-d41d-4c45-b7ca-d7f3cc168829
# ╠═c65e0994-8085-11ef-3cbc-fb515d8af323
# ╟─3325eff5-96ed-4823-b87e-185e75b3db97
# ╠═830393f6-65c8-49a8-bafc-447f142d04f3
# ╠═be615e70-389b-47a2-9a33-6b82dda1e650
# ╠═cdcd8c91-40d2-40f6-8865-306f45e6008f
# ╟─c3d7d282-9e6e-4c8f-885e-71bdebe210de
# ╠═c9ea9073-e973-42ba-b412-8a9938b11a25
# ╟─9ae7c1f5-a2bd-4518-b686-d53b81ad9fcb
# ╟─d45f9356-94f7-4320-ab75-24ef0b18fc80
# ╠═7469fefc-b084-4dc1-9c1a-78f747bd6eb0
# ╠═c4d133d7-5a9e-44b9-abcb-c544e0979135
# ╟─f54bbdb9-d9e2-47b0-aafe-5302b8510290
# ╟─1ca76166-5719-470c-8cd4-c44438310c86
# ╠═37e70e5a-eadd-4e3f-b6dd-e3ac05595f7e
# ╟─2925f271-54d7-4b83-a653-b157ac61f171
# ╟─04cc1bc5-2586-4102-914e-c1918b8868a7
# ╠═1c37cc21-b4b8-4ffb-92d9-696de4cfa462
# ╟─7df487d7-ca05-4fa8-a7ea-62298e25d400
# ╟─1f8b2a1f-4c51-46cf-9d0e-86bf1a6dba0b
# ╠═8a2f3599-c7be-4159-b56d-315b8f9b1996
# ╠═ff6559ef-9b86-40df-9ad6-4b921d934a6f
# ╟─4012a87c-03d6-45fe-b08e-41e0c3fb28a0
# ╠═27ce24a3-dca7-429e-849c-9f8accda58e4
# ╠═1ad3c07e-729c-480b-89b4-54e1f0f38f55
# ╟─42052bfb-8321-45c4-af8e-3d053e2ba623
# ╟─8d07f251-9362-495d-b502-c7b98a318140
# ╟─883fcfc6-bc40-4011-ada5-031259f42320
# ╠═37630f68-acd6-4517-a8c8-58f825b00c93
# ╠═d2da6ed8-d19c-4b12-ab33-384d0b524cb3
# ╠═fb0dbfdf-d2ce-4b20-bdeb-fa7cb6599378
# ╠═3b817455-d43d-4f61-ac3d-bdb80e992e34
# ╟─d6ba4c07-6c86-48de-979e-e212f22da557
# ╟─445da383-5dfe-427e-89f9-28ab480ffa47
# ╠═e566b981-c20d-4405-bbd3-7f59ecdcc915
# ╠═b3551f27-a745-4816-a173-6fcc54a2f001
# ╠═46f48e0e-e324-439c-9b71-762989b2b778
# ╟─87913ced-c8d0-42c1-916a-9250cdde8f73
# ╠═5eff0a2d-2e20-4699-a29f-d1be0cc4c428
# ╠═5846c57c-171b-4ef5-abcf-666596605e6d
# ╟─2d48d193-b822-46e2-a605-397b27634b3e
