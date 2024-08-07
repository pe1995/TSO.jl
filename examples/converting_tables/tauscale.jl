### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 19446b34-3546-11ef-3817-2354961c536a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using Plots
end

# ╔═╡ 24ef3a87-7f5d-4c8a-8eda-cfd9fb5b619d


# ╔═╡ 06cf8fae-b384-4938-9b33-ba4967b4abda
md"# Rosseland opacity tables"

# ╔═╡ 3050e202-654f-440d-96f6-41e0c23232d5
begin
    mother_table_path = "../../../opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"
    extension = "magg_m0_a0"
    eos_path = "ross_combined_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"
end;

# ╔═╡ 81162a6c-af24-4efd-b5a9-9fafc11c1dab
eos = reload(SqEoS, joinpath(mother_table_path, eos_path))

# ╔═╡ ef56d241-f292-4565-9fb2-554ed3049dbe
opa = reload(SqOpacity, joinpath(mother_table_path, opa_path), mmap=true)

# ╔═╡ 2c4db690-fb7d-456d-8775-8d30cd5f96e2
sun = @optical Average3D(eos, "../binning_opacities/sun_stagger.dat") eos opa

# ╔═╡ cb0915aa-d4cc-4055-81ff-12b107df3d14


# ╔═╡ bc9fb114-28c0-450c-afe4-98bf4aac3b57
md"# Interpolation to reference opacity"

# ╔═╡ 7fb1eb25-412e-4d0c-8bd8-8b595cf8ff85
λ_ref = 5000.0

# ╔═╡ d9a4d561-6904-4fa0-bda5-0a1c82a1bee4
begin
	κ_ref = zeros(size(eos))
	κ_i = zeros(length(opa.λ))
	for j in axes(opa.κ, 2)
		for i in axes(opa.κ, 1)
			κ_i .= log.(opa.κ[i, j, :])
			κ_ref[i, j] = TSO.linear_interpolation(opa.λ, κ_i)(λ_ref)
		end
	end
end

# ╔═╡ 1bb6829d-c3e8-4a16-9997-3d222b6e0c8a
eos2 = deepcopy(eos)

# ╔═╡ 3e33a373-962b-46bd-868f-4fc4f65a2eb9
eos2.lnRoss .= κ_ref

# ╔═╡ 557c5b97-f300-4069-b6ad-a7407299dc19


# ╔═╡ 5a7a72f4-b389-4c25-9756-12cae9b6eb23
sun_ref = @optical Average3D(eos2, "../binning_opacities/sun_stagger.dat") eos2

# ╔═╡ d571492c-0ecf-4a46-a4bb-ca2cab15c58f
let
	plot(framestyle=:box, grid=false)
	plot!(log10.(sun.τ), exp.(sun.lnT), label="Rosseland")
	plot!(log10.(sun_ref.τ), exp.(sun_ref.lnT), label="τ500")
end

# ╔═╡ 48e10d6b-763e-4037-931e-869f70875356


# ╔═╡ 4e1ddfd3-5cc1-4c8f-8663-e00d9118d286
save(eos2, joinpath(mother_table_path, "tau500_combined_eos_"*extension*".hdf5"))

# ╔═╡ Cell order:
# ╠═19446b34-3546-11ef-3817-2354961c536a
# ╟─24ef3a87-7f5d-4c8a-8eda-cfd9fb5b619d
# ╟─06cf8fae-b384-4938-9b33-ba4967b4abda
# ╠═3050e202-654f-440d-96f6-41e0c23232d5
# ╠═81162a6c-af24-4efd-b5a9-9fafc11c1dab
# ╠═ef56d241-f292-4565-9fb2-554ed3049dbe
# ╠═2c4db690-fb7d-456d-8775-8d30cd5f96e2
# ╟─cb0915aa-d4cc-4055-81ff-12b107df3d14
# ╟─bc9fb114-28c0-450c-afe4-98bf4aac3b57
# ╠═7fb1eb25-412e-4d0c-8bd8-8b595cf8ff85
# ╠═d9a4d561-6904-4fa0-bda5-0a1c82a1bee4
# ╠═1bb6829d-c3e8-4a16-9997-3d222b6e0c8a
# ╠═3e33a373-962b-46bd-868f-4fc4f65a2eb9
# ╟─557c5b97-f300-4069-b6ad-a7407299dc19
# ╠═5a7a72f4-b389-4c25-9756-12cae9b6eb23
# ╠═d571492c-0ecf-4a46-a4bb-ca2cab15c58f
# ╟─48e10d6b-763e-4037-931e-869f70875356
# ╠═4e1ddfd3-5cc1-4c8f-8663-e00d9118d286
