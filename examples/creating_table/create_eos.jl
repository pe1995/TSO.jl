### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1ebcbc84-20c4-11ee-2210-8b27cd286dfb
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using TSO
	using Plots
	
	Glob = TSO.Glob;
end

# ╔═╡ 4de1c8c8-e44b-4305-bafb-0a34d0d911e1
md"# Turbospectrum Equation of State
The EoS that is used in DISPATCH can be computed from first principles using the Turbospectrum code. Turbospectrum can be used to compute opacities as well, although this is currently not the main channel. At the moment, the code only uses the EoS from this script, and complements it with opacities from the MARCS code."

# ╔═╡ 6284d0e6-917a-49a2-ac0f-6912b1d67408
md"For the setup we use the Turspectrum Wrapper, although for the actual execution the wrapper is not used, it is a julia-specific wrapper"

# ╔═╡ 60c90eec-8f36-47bd-9320-e198dab6fd70
begin
	TSO.load_TS()         # Set the TS path
	TSO.load_wrapper()    # Set the Wrapper path
	TSO.import_wrapper()  # import the python modules of the wrapper
end

# ╔═╡ 12b4c825-766b-4cbd-acc6-b5998c6247c6
!isdir(TSO.@inWrapper("example/models")) && mkdir(TSO.@inWrapper("example/models"))

# ╔═╡ 9ae940b1-27f4-46e0-9707-299df4b9b7b7
md"## The Grid"

# ╔═╡ a128a520-895b-4f40-9fbb-6586ac25da10
lnT = range(log(1.1e3), log(5.5e5); length=159)

# ╔═╡ af66c1ea-fee0-4172-9100-88ad54993ebd
lnρ = range(log(1e-15), log(1e-3); length=159)

# ╔═╡ 5a695141-8e7c-4b85-a886-9c42e0fdad55
TSO.write_as_stagger(Float64[lnT...], Float64[lnρ...])

# ╔═╡ a420dd22-8e80-451d-9fe8-f24f161bdfb6
md"## The Setup"

# ╔═╡ 5b88ccf9-928d-4ecc-898e-cc15122ecbe3
begin
	## config file used to setup TurboSpectrum (TS) run
    debug = 1

    ## TS root directory
    ts_root = TSO.@inTS ""
   
    ## list of model atmospheres to use
    atmos_list = TSO.@inWrapper "example/models/TSO_list.in"

    ## path to model atmospheres 
    atmos_path = TSO.@inWrapper "example/models/"

    ## 'm1d' or 'marcs' -> always Stagger for opacity tables
    atmos_format = "stagger"
end

# ╔═╡ 76e84b2d-b6fc-473a-8a32-850c2d23ae7c
linelist = abspath.([
		"LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list", 
        "LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
        "LINE-LISTS/ADDITIONAL-LISTS/Hlinedata"
	]
)

# ╔═╡ 95027847-2d31-43e9-84e2-576e42e38c36
begin
	lam_start = 1000
	lam_end   = 4000
	
	## resolution per wavelenght (R capital)
	resolution = 200000 # hres
	#resolution = 20000  # lres
end

# ╔═╡ 6a727aa9-1779-4d5d-a25d-dcf2d52d17c6
tmolim = 20000.0

# ╔═╡ 13023f21-8428-4519-b0d1-59adc2d989ad
@info "Chosen λ step + Number of points: $(TSO.ΔΛ(lam_start,lam_end,resolution)), $(TSO.N_Λ(lam_start,lam_end,resolution))"

# ╔═╡ 7a5f3c16-095b-4ef1-9832-30b24066d4da
setup_input = Dict(
	"debug"        =>debug, 
	"ts_root"      =>ts_root, 
	"atmos_path"   =>atmos_path, 
	"atmos_format" =>atmos_format, 
	"atmos_list"   =>atmos_list, 
	"linelist"     =>TSO.Py(String[]).to_numpy(),  
	"lam_start"    =>lam_start, 
	"lam_end"      =>lam_end, 
	"resolution"   =>resolution,
	"TMOLIM"       =>tmolim
)

# ╔═╡ 7ac5674c-b5d9-43dd-a0c7-64b269e26295
begin
	## Create the setup object
	setup       = TSO.computeOpac.setup(file=setup_input, mode="MAprovided")
	setup.jobID = "TSO"
	
	wvl_set = "asplund07_v5.0"
end

# ╔═╡ 04c77365-c905-4635-b9a2-dc5e44628a23
begin
    magg_2022 = [
		(id, TSO.magg2022_abund[String(TSO.id_atom(id))]) 
						for id in eachindex(TSO.atomic_number)
							if String(TSO.id_atom(id)) in keys(TSO.magg2022_abund) 
	]

	asplund_2007 = [
		(id, TSO.asplund2007_abund[String(TSO.id_atom(id))]) 
						for id in eachindex(TSO.atomic_number)
							if String(TSO.id_atom(id)) in keys(TSO.asplund2007_abund) 
	]
		
    # Pick the abundances
    abundances = asplund_2007
    @info "Modify the following abundances: (Species, ID, Abundance)"
    for a in abundances
        @info "$(TSO.id_atom(a[1]))-$(a[1]) -> $(a[2])"
    end
end

# ╔═╡ 3ec81df5-234b-4315-8156-00a5f8d8c5fb
md"## Running Turbospectrum
Turbospectrum can be run in 2 steps, for the EoS only the quick first step is important. We can run this step easily within an slurm allocation as job steps. The abundances can be chosen in the input setup from the wrapper. If the following line is executed without a slurm allocation, there is no memory management at the moment!"

# ╔═╡ 8e6783e6-259e-4829-bd34-18d1d9b86268
TSO.babsma!(setup, abundances)

# ╔═╡ 59861ad5-27ec-47b5-9d18-fd812cb1b58c
md"## Collect the output"

# ╔═╡ 2021e20f-6699-45e7-951e-4795a70b206c
lot = Glob.glob("_TSOeos_*_TSO.eos")

# ╔═╡ 7e497fcd-2a3e-4d5e-be8a-f8060bf7f131
eos = TSO.load(TSO.EoSTable, lot) 

# ╔═╡ e39c4a02-266a-4038-8de1-9604f1c967d8
md"Apply smoothing if needed"

# ╔═╡ 54e669b2-bb26-450f-a0e8-b2fd1fc7706d
TSO.smooth!(eos)

# ╔═╡ 338f7790-a7e5-4961-9230-56e66debdbcb
begin
	xx, yy = meshgrid(eos.lnT, eos.lnRho)
	scatter(xx, yy, marker_z=eos.lnEi, 
		markersize=2, 
		label=nothing, 
		colormap=:rainbow, markerstrokewidth=0)
end

# ╔═╡ cf3a11ab-cb6d-4ecc-8346-5a44c3986eeb
save(eos, "eos_$(wvl_set).hdf5")

# ╔═╡ 4b0defa3-e561-4627-8b98-f1c2f9b916bf
TSO.move_output() # clean up

# ╔═╡ Cell order:
# ╟─4de1c8c8-e44b-4305-bafb-0a34d0d911e1
# ╠═1ebcbc84-20c4-11ee-2210-8b27cd286dfb
# ╟─6284d0e6-917a-49a2-ac0f-6912b1d67408
# ╠═60c90eec-8f36-47bd-9320-e198dab6fd70
# ╠═12b4c825-766b-4cbd-acc6-b5998c6247c6
# ╟─9ae940b1-27f4-46e0-9707-299df4b9b7b7
# ╠═a128a520-895b-4f40-9fbb-6586ac25da10
# ╠═af66c1ea-fee0-4172-9100-88ad54993ebd
# ╠═5a695141-8e7c-4b85-a886-9c42e0fdad55
# ╟─a420dd22-8e80-451d-9fe8-f24f161bdfb6
# ╠═5b88ccf9-928d-4ecc-898e-cc15122ecbe3
# ╠═76e84b2d-b6fc-473a-8a32-850c2d23ae7c
# ╠═95027847-2d31-43e9-84e2-576e42e38c36
# ╠═6a727aa9-1779-4d5d-a25d-dcf2d52d17c6
# ╠═13023f21-8428-4519-b0d1-59adc2d989ad
# ╠═7a5f3c16-095b-4ef1-9832-30b24066d4da
# ╠═7ac5674c-b5d9-43dd-a0c7-64b269e26295
# ╠═04c77365-c905-4635-b9a2-dc5e44628a23
# ╟─3ec81df5-234b-4315-8156-00a5f8d8c5fb
# ╠═8e6783e6-259e-4829-bd34-18d1d9b86268
# ╟─59861ad5-27ec-47b5-9d18-fd812cb1b58c
# ╠═2021e20f-6699-45e7-951e-4795a70b206c
# ╠═7e497fcd-2a3e-4d5e-be8a-f8060bf7f131
# ╟─e39c4a02-266a-4038-8de1-9604f1c967d8
# ╠═54e669b2-bb26-450f-a0e8-b2fd1fc7706d
# ╠═338f7790-a7e5-4961-9230-56e66debdbcb
# ╠═cf3a11ab-cb6d-4ecc-8346-5a44c3986eeb
# ╠═4b0defa3-e561-4627-8b98-f1c2f9b916bf
