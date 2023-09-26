### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 7cec369a-e99a-11ed-0d22-53109bf030b5
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
	#using PyPlot
	#using Dierckx
	using Glob
	using Serialization
end

# ╔═╡ 022571bd-bdd4-4096-a94e-03f59ea850f7
md"
# Reading Raw Tables
Pick the EoS that shall provide the internal energy to the opacity table.
"

# ╔═╡ f02fbdc0-1d14-4c5e-8ce8-091b07d083da
aos = @axed reload(SqEoS, abspath("../creating_table/eos_asplund07_m1_a04_v5.0.hdf5"))

# ╔═╡ c7f9ebfc-4010-47e5-afe0-183efec4273f
md"Folder of the Opacity tables"

# ╔═╡ d5be07a6-c7bc-48ea-8337-73f445dfea1a
paths = glob("OS_table*", "opacity_tables/MARCS/asplund/Z-1.0a0.4/")

# ╔═╡ 626c6ea3-ee30-4dbe-9e32-4b211fa559e0
md"Read the raw opacity tables"

# ╔═╡ 75d03d5a-ef5f-47c6-8a5c-f14ce7cba955
mos = MARCSOpacity(paths...)

# ╔═╡ d8634756-2bd8-4e9f-b8d0-2d45ed301495
md"The raw tables look like the following"

# ╔═╡ f63561ff-b4e9-4ddc-b5f5-f69a24dc11b9
begin
	#=close()

	for i in eachindex(mos)
		k = [minimum(mos[i].κ_c[j, :]) for j in axes(mos[i].κ_c, 1)]
	    plt.scatter(log.(mos[i].T), log.(mos[i].ρ), c=log10.(mos[i].κ_c[:, 1000]), 
			s=3, vmin=0.5, vmax=3.5)
	end
	
	c = plt.colorbar()
	c.set_label("minimum opacity")
	gcf()=#
end

# ╔═╡ 03a23793-7af7-4e6b-9a20-33452dc5f6a8
for i in eachindex(mos)
	@show minimum(mos[i].κ_c[:, 1])
end

# ╔═╡ e379074f-997a-45a3-922d-8e57b7709502
md"# Experiment: Unstructured Interpolation"

# ╔═╡ 1a353858-af6b-447b-8023-045f73c47b46
#m_int = uniform(mos..., new_T_size=104, new_ρ_size=104, structured=false)

# ╔═╡ 276ea20d-ab1e-44ca-bcfc-59c34f936bdb
md"
# Interpolate
We can now first interpolate the unstructured data to a square T-rho table, which
is required for the conversion to dispatch later, and also to fit into the rest of the API."

# ╔═╡ 18876950-d744-40c7-9b8c-1ec8044708f9
m_int = uniform(mos..., new_T_size=104, new_ρ_size=104)

# ╔═╡ 772769bf-d366-4265-aa72-91ba2ee5d803
begin
	#=close()


	tt, rr = meshgrid(m_int.T, m_int.ρ)
	
	plt.scatter(tt, rr, c=log10.(m_int.κ_c[:, :, 1000]), s=2, vmin=0.5, vmax=3.5)
	
	c2 = plt.colorbar()
	c2.set_label("minimum line opacity")

	
	gcf()=#
end

# ╔═╡ c6585349-3970-440d-ab58-45d82ef5943e
@show minimum(m_int.κ_c[:, :, 1])

# ╔═╡ a1457691-0d92-45b9-b625-49dfc3a349c2
md"Next, we add the internal energy from an external source. In this case, we add it from an independent Turbospectrum run."

# ╔═╡ 2cf4938f-7d0e-474e-90cf-e6c7d13e3e36
neweos, newopa, newopa_c, newopa_l, newopa_s = complement(m_int, aos, unify=false)

# ╔═╡ 486d3b00-b2ca-450e-ba01-6ceeac3c3ac8
md"We set the limits of the tables to be 1e±30"

# ╔═╡ 2696a129-2c29-4e06-93bf-a8240a76d6c9
begin
	set_limits!(@axed(neweos), newopa)
	set_limits!(@axed(neweos), newopa_c)
	set_limits!(@axed(neweos), newopa_l)
end

# ╔═╡ 36d6edb6-0ec8-4169-bd17-44bb66fea284
md"
# Save and Move
Save everything in the usual TSO.jl format
"

# ╔═╡ fd21fd34-a8a1-4f14-8e54-a3cbe76a17cb
begin
	dname = "TSO_MARCS_asplund_m1_a04_v1.6"

	if !isdir(dname)
		mkdir(dname)
	else
		rm(dname, recursive=true)
		mkdir(dname)
	end
	
	save(neweos,   joinpath(dname, "combined_eos.hdf5"))
	save(newopa,   joinpath(dname, "combined_opacities.hdf5"))
	save(newopa_c, joinpath(dname, "combined_Copacities.hdf5"))
	save(newopa_l, joinpath(dname, "combined_Lopacities.hdf5"))
	save(newopa_s, joinpath(dname, "combined_Sopacities.hdf5"))
end

# ╔═╡ Cell order:
# ╠═7cec369a-e99a-11ed-0d22-53109bf030b5
# ╟─022571bd-bdd4-4096-a94e-03f59ea850f7
# ╠═f02fbdc0-1d14-4c5e-8ce8-091b07d083da
# ╟─c7f9ebfc-4010-47e5-afe0-183efec4273f
# ╠═d5be07a6-c7bc-48ea-8337-73f445dfea1a
# ╟─626c6ea3-ee30-4dbe-9e32-4b211fa559e0
# ╠═75d03d5a-ef5f-47c6-8a5c-f14ce7cba955
# ╟─d8634756-2bd8-4e9f-b8d0-2d45ed301495
# ╠═f63561ff-b4e9-4ddc-b5f5-f69a24dc11b9
# ╠═03a23793-7af7-4e6b-9a20-33452dc5f6a8
# ╟─e379074f-997a-45a3-922d-8e57b7709502
# ╠═1a353858-af6b-447b-8023-045f73c47b46
# ╟─276ea20d-ab1e-44ca-bcfc-59c34f936bdb
# ╠═18876950-d744-40c7-9b8c-1ec8044708f9
# ╟─772769bf-d366-4265-aa72-91ba2ee5d803
# ╠═c6585349-3970-440d-ab58-45d82ef5943e
# ╟─a1457691-0d92-45b9-b625-49dfc3a349c2
# ╠═2cf4938f-7d0e-474e-90cf-e6c7d13e3e36
# ╟─486d3b00-b2ca-450e-ba01-6ceeac3c3ac8
# ╠═2696a129-2c29-4e06-93bf-a8240a76d6c9
# ╟─36d6edb6-0ec8-4169-bd17-44bb66fea284
# ╠═fd21fd34-a8a1-4f14-8e54-a3cbe76a17cb
