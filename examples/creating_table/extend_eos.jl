### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 14822926-c01e-11ee-38e1-13e9ffda71eb
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end

# ╔═╡ f100dee1-f01d-4c15-bd7c-ac56afc76b9a


# ╔═╡ 7b912e88-d9fb-4b7f-b97d-5c339b45adfe
md"Read the default tables"

# ╔═╡ 4c84a8ad-168e-4386-a993-644cae35843f
begin
	oppath = "../../../MUST.jl/examples/initial_models/DIS_MARCS_t50g40m00_v0.5"
	eos = reload(SqEoS, joinpath(oppath, "eos.hdf5"))
	opa = reload(SqOpacity, joinpath(oppath, "binned_opacities.hdf5"))
end

# ╔═╡ 255201cb-2d2d-4b7b-bc14-1453ceddeb3e


# ╔═╡ 56eec877-d6b1-4c45-b5c2-d46a320f55b6
md"We can extend those tables directly using the `extend()` function"

# ╔═╡ f44b14d9-e12e-46e6-b2a0-32fc046625d2
eos_new, opa_new = TSO.extend(eos, opa, downD=0.2, downE=0.2, upD=0.1, upE=0.1)

# ╔═╡ d30e31e8-7aea-4720-9e04-fb3d0182a37f
begin
	TSO.set_limits!(@axed(eos_new), opa_new)
	TSO.fill_nan!(@axed(eos_new), opa_new)
end

# ╔═╡ e6f76b32-1c10-4fab-893b-c056cc7c0f71


# ╔═╡ 5677f4d0-cf78-4b8a-9637-3b06190ff7ab
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	xx, yy = meshgrid(@axed(eos))
	
	im = ax.scatter(xx, yy, s=0.1, c=log10.(opa.κ[:, :, 1]))
	plt.colorbar(im, ax=ax)

	ax.set_xlim(limits(eos_new)[1:2]...)
	ax.set_ylim(limits(eos_new)[3:4]...)
	
	gcf()
end

# ╔═╡ aabfd63f-f115-4bc0-8325-372e3664ebc4
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	xx, yy = meshgrid(@axed(eos_new))
	
	im = ax.scatter(xx, yy, s=0.1, c=log10.(opa_new.κ[:, :, 1]))
	plt.colorbar(im, ax=ax)

	ax.set_xlim(limits(eos_new)[1:2]...)
	ax.set_ylim(limits(eos_new)[3:4]...)
	
	gcf()
end

# ╔═╡ bb820023-bd9f-4729-a021-da7f857b74c4


# ╔═╡ a3a912e0-dddb-4606-b06b-1a5a0b2c3d92
@info "old size <-> new size: $(size(eos)) <-> $(size(eos_new))"

# ╔═╡ Cell order:
# ╠═14822926-c01e-11ee-38e1-13e9ffda71eb
# ╟─f100dee1-f01d-4c15-bd7c-ac56afc76b9a
# ╟─7b912e88-d9fb-4b7f-b97d-5c339b45adfe
# ╠═4c84a8ad-168e-4386-a993-644cae35843f
# ╟─255201cb-2d2d-4b7b-bc14-1453ceddeb3e
# ╟─56eec877-d6b1-4c45-b5c2-d46a320f55b6
# ╠═f44b14d9-e12e-46e6-b2a0-32fc046625d2
# ╠═d30e31e8-7aea-4720-9e04-fb3d0182a37f
# ╟─e6f76b32-1c10-4fab-893b-c056cc7c0f71
# ╟─5677f4d0-cf78-4b8a-9637-3b06190ff7ab
# ╟─aabfd63f-f115-4bc0-8325-372e3664ebc4
# ╟─bb820023-bd9f-4729-a021-da7f857b74c4
# ╟─a3a912e0-dddb-4606-b06b-1a5a0b2c3d92
