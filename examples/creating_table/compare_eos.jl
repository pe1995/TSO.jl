### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 9ec4d44c-2165-11ee-0787-693b454dec92
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
	using Plots
end

# ╔═╡ 0458da53-5933-46de-8deb-917b3daaa868
eos_magg = reload(SqEoS, "eos_magg_v5.0.hdf5")

# ╔═╡ 52723abb-db1f-4883-856d-3b917497bc5d
eos_asplund = reload(SqEoS, "eos_asplund07_v5.0.hdf5")

# ╔═╡ ad03031a-92b2-4c79-9ece-70e7d27a88d4
begin
	xx, yy = meshgrid(eos_magg.lnT, eos_magg.lnRho)
	
	scatter(xx, yy, marker_z=100 .*(eos_magg.lnEi.-eos_asplund.lnEi) ./ eos_asplund.lnEi, 
		markersize=2, 
		label=nothing, 
		colormap=:rainbow, markerstrokewidth=0)
end

# ╔═╡ Cell order:
# ╠═9ec4d44c-2165-11ee-0787-693b454dec92
# ╠═0458da53-5933-46de-8deb-917b3daaa868
# ╠═52723abb-db1f-4883-856d-3b917497bc5d
# ╠═ad03031a-92b2-4c79-9ece-70e7d27a88d4
