### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 41cfb69a-2173-11ee-047c-4ff531cb20fe
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
end

# ╔═╡ 64291703-4643-4856-b016-ee8e5f06759f
md"# Complement a given opacity table with a give EoS"

# ╔═╡ 5eeb2e09-4741-4598-ae52-a5172db25d06
extension = "_magg22"

# ╔═╡ 521ea7cb-cd8f-43cc-a1ca-b8081a4db097
aos = @axed reload(SqEoS, abspath("../creating_table/eos_magg_v5.0.hdf5"))

# ╔═╡ d0b263fd-7709-4499-b139-988022ebb851
folderopa = abspath("../../../opacity_tables/TSO_MARCS_v1.6")

# ╔═╡ 69024906-96c2-4548-9a16-1ba561bfa2ee
eos_old = reload(SqEoS, joinpath(folderopa, "combined_eos.hdf5"), mmap=true)

# ╔═╡ 1e79c58b-2c02-4c8a-9aee-ee68ea2d3d97
opa = reload(SqOpacity, joinpath(folderopa, "combined_opacities.hdf5"), mmap=true)

# ╔═╡ ef60e67f-9204-421b-8a47-0fb1d124bca9
opas = reload(SqOpacity, joinpath(folderopa, "combined_Sopacities.hdf5"), mmap=true)

# ╔═╡ e69a9bce-7317-4ece-8c6a-0be2689789ca
opal = reload(SqOpacity, joinpath(folderopa, "combined_Copacities.hdf5"), mmap=true)

# ╔═╡ 66ffb5fb-1c4b-4323-84d2-907ac6a4d408
opac = reload(SqOpacity, joinpath(folderopa, "combined_Lopacities.hdf5"), mmap=true)

# ╔═╡ 74bc181c-51c3-4a6c-b313-0b6cbf5a26d0
md"## Interpolate opacities to the grid of the EoS"

# ╔═╡ 3c50a4d4-d4bd-4dc0-86fb-bd923fb4f7b3
md"Are they both on the same energy scale?"

# ╔═╡ 845b4944-8866-4e27-9910-e5c21146b2b8
is_internal_energy(@axed eos_old) == is_internal_energy(aos)

# ╔═╡ 449bd01b-de51-44f7-99f9-a9eb6be1d938
begin
	opa_new, opas_new = complement(@axed(eos_old), aos, opa, opas)
	#opas_new = complement(@axed(eos_old), aos, opas)

	set_limits!(aos, opa_new)
	set_limits!(aos, opas_new)
end

# ╔═╡ 72424b05-8d9d-46d3-9f54-ed4b45e624b4
begin
	save(opa_new, joinpath(folderopa, "combined_opacities$(extension).hdf5"))
	save(opas_new, joinpath(folderopa, "combined_Sopacities$(extension).hdf5"))
end

# ╔═╡ f4e7dae8-c8c6-4963-a41a-aac6d633463e
save(aos.eos, joinpath(folderopa, "combined_eos$(extension).hdf5"))

# ╔═╡ Cell order:
# ╟─64291703-4643-4856-b016-ee8e5f06759f
# ╠═41cfb69a-2173-11ee-047c-4ff531cb20fe
# ╠═5eeb2e09-4741-4598-ae52-a5172db25d06
# ╠═521ea7cb-cd8f-43cc-a1ca-b8081a4db097
# ╠═d0b263fd-7709-4499-b139-988022ebb851
# ╠═69024906-96c2-4548-9a16-1ba561bfa2ee
# ╠═1e79c58b-2c02-4c8a-9aee-ee68ea2d3d97
# ╠═ef60e67f-9204-421b-8a47-0fb1d124bca9
# ╠═e69a9bce-7317-4ece-8c6a-0be2689789ca
# ╠═66ffb5fb-1c4b-4323-84d2-907ac6a4d408
# ╟─74bc181c-51c3-4a6c-b313-0b6cbf5a26d0
# ╟─3c50a4d4-d4bd-4dc0-86fb-bd923fb4f7b3
# ╠═845b4944-8866-4e27-9910-e5c21146b2b8
# ╠═449bd01b-de51-44f7-99f9-a9eb6be1d938
# ╠═72424b05-8d9d-46d3-9f54-ed4b45e624b4
# ╠═f4e7dae8-c8c6-4963-a41a-aac6d633463e
