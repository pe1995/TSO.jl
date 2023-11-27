### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 3c3d9886-587b-4209-a421-1f21e6749e04
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using DelimitedFiles
end

# ╔═╡ 317e3358-8512-4166-9596-23ae7bb380c6
md"# Converting raw opacity tables to DISPATCH"

# ╔═╡ 48efb4ac-ca3f-4a48-b5ed-9e173fb01846
md"## Functionality"

# ╔═╡ 2de77db3-f110-4e20-b68e-b00ce9451e46
"""
	formation_opacity(table_folder, av_path; logg=log10(2.75e4))

Compute the rosseland optical depth for the entire table. Afterwards, compute the formation height of every wavelength point in the unbinned opacity table.
"""
function formation_opacities(table_folder, av_path; logg=log10(2.75e4), extension="")

    opacities = reload(SqOpacity, 
				joinpath(table_folder, "combined_opacities$(extension).hdf5"))
    eos = reload(SqEoS,     
				joinpath(table_folder, "combined_eos$(extension).hdf5"))
    aos = @axed eos


    model = Average3D(av_path, logg=logg)

  
    if !isfile(joinpath(table_folder, "combined_ross_eos$(extension).hdf5"))
        @info "computing rosseland"
        rosseland_opacity!(aos.eos.lnRoss, aos, opacities; 
					weights=ω_midpoint(opacities))
		
        transfer_rosseland!(aos.eos, opacities)
		
        save(aos.eos,   joinpath(table_folder, "combined_ross_eos$(extension).hdf5"))
        save(opacities, joinpath(table_folder, "combined_opacities$(extension).hdf5"))
    end


    if isfile(joinpath(table_folder, "combined_formation_opacities$(extension).hdf5"))
        @warn "skipping formation opacities"
        return
    end

    τ_ross, τ_λ = optical_depth(aos, opacities, model)
    d_ross, d_κ = formation_height(model, aos, opacities, τ_ross, τ_λ)


    formation_opacities = SqOpacity(d_κ, d_ross, opacities.src, opacities.λ, true)
	
    save(formation_opacities, 
			joinpath(table_folder, "combined_formation_opacities$(extension).hdf5"))

end

# ╔═╡ aad8b35a-ea58-11ed-2350-79fec9c7f69c
"""
	bin_opacities(table_folder, av_path, new_folder; logg=log10(2.75e4))

Bin opacities, based on a specific binning method, using formation opacities computed previously. Save the binned tables.
"""
function bin_opacities(table_folder, av_path, name; 
                method=:uniform, 
                use_contribution=false, 
                logg=log10(2.75e4), 
                Nbins=5, 
                line_bins=3, 
                extension="",
                kwargs...)

    eos = reload(SqEoS, 
			joinpath(table_folder, "combined_ross_eos$(extension).hdf5")) 
	
    opacities = reload(SqOpacity, 
			joinpath(table_folder, "combined_opacities$(extension).hdf5"), 
			mmap=true)

    sopacities = reload(SqOpacity, 
			joinpath(table_folder, "combined_Sopacities$(extension).hdf5"), 
			mmap=true)
	
    formOpacities = reload(SqOpacity, 
			joinpath(table_folder, "combined_formation_opacities$(extension).hdf5"), 
			mmap=true)

    #κ_s = sopacities.κ
    #κ = opacities.κ
    #@show minimum(κ_s) maximum(κ_s) any(isnan.(κ_s)) all(isfinite.(κ_s))
    #@show minimum(κ) maximum(κ) any(isnan.(κ)) all(isfinite.(κ))


    weights = ω_midpoint(opacities)
    model = @optical Average3D(eos, av_path, logg=logg) eos opacities


    TSO.Clustering.Random.seed!(42)


    ks = if use_contribution
        T = TSO.linear_interpolation(log10.(model.τ), model.lnT, extrapolation_bc=TSO.Line()).(log10.(formOpacities.κ_ross))
        S = TSO.Bλ.(opacities.λ, Base.convert.(Float32, exp.(T)))
        -log10.(formOpacities.κ_ross .* S)
    else
        -log10.(formOpacities.κ_ross)
    end


    bins = if method == :tabgen
        TabgenBinning(TSO.UniformTabgenBins; 
                    opacities=opacities, 
                    formation_opacity=ks, 
                    Nbins=Nbins, line_bins=line_bins, kwargs...)

    elseif method == :equal
        TabgenBinning(TSO.EqualTabgenBins;
                    opacities=opacities, 
                    formation_opacity=ks, 
                    Nbins=Nbins, kwargs...)
    elseif method == :stagger
        StaggerBinning(TSO.StaggerBins;
                    opacities=opacities, 
                    formation_opacity=ks, 
                    Nbins=Nbins, kwargs...)
    elseif method == :kmeans
        ClusterBinning(TSO.KmeansBins;
                            opacities=opacities, 
                            formation_opacity=ks, 
                            Nbins=Nbins, kwargs...)
    else
        error("Chosen method not implemented.")
    end

    bin = binning(bins, opacities, ks) 





    function save_table(binned_opacities)
        eos_table_name = name

        !isdir(eos_table_name) && mkdir(eos_table_name) 

        cp(joinpath(table_folder, "combined_ross_eos$(extension).hdf5"), 
			joinpath(eos_table_name, "eos.hdf5"), force=true)

        save(binned_opacities.opacities, joinpath(eos_table_name, "binned_opacities.hdf5"))

        fid = TSO.HDF5.h5open(joinpath(eos_table_name, "bin_assignment.hdf5"), "w")
        fid["bins"] = bin
        fid["lambda"] = opacities.λ
        close(fid)
    end
	
    binned_opacities = tabulate(bin, weights, eos, opacities, sopacities, transition_model=model)
    #binned_opacities = tabulate(bin, weights, eos, opacities, transition_model=model)
    
    save_table(binned_opacities)
end

# ╔═╡ 9aae590d-b656-43f3-949e-0966f1a7a32a
"""
	fromT_toE(table_folder, folder_new)

Convert the binned opacities + eos from a T-ρ to a Eint-ρ grid. Upsample the resulting table.
"""
function fromT_toE(table_folder, folder_new; upsample=1000)
    eos = reload(SqEoS,     joinpath(table_folder, "eos.hdf5"))
    opa = reload(SqOpacity, joinpath(table_folder, "binned_opacities.hdf5"));

    aos = @axed eos
    eosE, opaE = switch_energy(aos, opa, upsample=upsample);
    aosE = @axed eosE

    TSO.fill_nan!(aosE, opaE)
    TSO.set_limits!(aosE, opaE)

    !isdir(folder_new) && mkdir(folder_new) 

    for_dispatch(eosE, opaE.κ, opaE.src, ones(eltype(opaE.src), size(opaE.src)...), folder_new)

    save(opaE, joinpath(folder_new, "binned_opacities.hdf5"))
    save(eosE, joinpath(folder_new, "eos.hdf5"))
end

# ╔═╡ 9fb3f203-331e-4b76-b94c-c992f625c80a
md"## Input parameters"

# ╔═╡ b45ecd12-4a70-42ac-bbee-9d6568adc214
begin
	eos_folder = "tables/TSO_MARCS_v1.6"
	model = "sun_stagger.dat"
	new_eos_folder = "DIS_MARCS_v1.7.4"
	new_eosE_folder = "DIS_MARCS_E_v1.7.4"
    extension = "_magg22"
end;

# ╔═╡ 0162409f-700e-457f-bd2d-ca685a025274
md"## Converting"

# ╔═╡ 94c6f318-22be-47e9-9248-adf5877e4d8a
formation_opacities(eos_folder, model, extension=extension)

# ╔═╡ a8c93f35-2ee0-4e62-b9f2-4762d108f02e
bin_opacities(eos_folder, model, new_eos_folder, 
        extension=extension,
        method=:kmeans, 
        use_contribution=false, 
        stripes=false,
        Nbins=16, 
        quadrants=[ 
            TSO.Quadrant((0.0, 4.0),   (-100, 4.5), 4),
            TSO.Quadrant((0.0, 4.0),   (4.5, 100), 4),
            TSO.Quadrant((4.0, 100.0), (-100, 100), 8)
        ],
            #=
            TSO.Quadrant((0.0, 3.6),  (-100, 0.5), 3), 
            TSO.Quadrant((0.0, 3.6),  (0.5, 4.0),  3), 
            TSO.Quadrant((0.0, 3.6),  (4.0, 100),  3), 

            TSO.Quadrant((3.6, 4.0),  (-100, 0.5), 3), 
            TSO.Quadrant((3.6, 4.0),  (0.5, 100), 3), 

            TSO.Quadrant((4.0, 4.5),  (-100, 0.5), 3), 
            TSO.Quadrant((4.0, 4.5),  (0.5, 100), 3),  

            TSO.Quadrant((4.5, 5.0),  (-100, 100), 3), 

            TSO.Quadrant((5.0, 100),  (-100, 100), 3)],
            =#

            #TSO.Quadrant((0.0, 3.55), (-100, 4.0), 5), 
            #TSO.Quadrant((0.0, 3.55), (4.0,  100), 1), 
            #TSO.Quadrant((3.55, 4.0), (-100, 100), 5), 
            ##TSO.Quadrant((3.7, 4.0),  (-100, 100), 5), 
            #TSO.Quadrant((4.0, 4.5),  (-100, 100), 5), 
            #TSO.Quadrant((4.5, 5.0),  (-100, 100), 1), 
            #TSO.Quadrant((5.0, 100),  (-100, 100), 1)],

            #=
            #TSO.Quadrant((0.0, 3.85),  (-100, 1.0), 3),   
            TSO.Quadrant((0.0, 3.55),  (-100, 5.0), 5),      # very important!
            TSO.Quadrant((0.0, 3.55),  (5.0, 100),  2),      # quite important far out
            #TSO.Quadrant((3.55, 3.80),  (-100, 1.0),  2),    # not as important as I thought
            #TSO.Quadrant((3.55, 3.80),  (1.0,  100),  3),    # quite important far out
            #TSO.Quadrant((3.55, 4.05),  (-100, 100.0),  5),  # very important!
            TSO.Quadrant((3.55, 3.70),  (-100, 100.0),  5),
            TSO.Quadrant((3.70, 4.05),  (-100, 100.0),  5),
            TSO.Quadrant((4.05, 4.55),  (-100, 100.0),  3),  # quite important
            #TSO.Quadrant((4.30, 4.55),  (-100, 100.0),  1),
            TSO.Quadrant((4.55, 4.80),  (-100, 100.0),  1),
            TSO.Quadrant((4.80, 100),   (-100, 100.0),  1),
            ##TSO.Quadrant((4.55, 100),  (1.0,  100),  1)
		],
        #=quadrants=[ 
            TSO.Quadrant((0.0, 3.6),   (-100, 4.5), 3),
            TSO.Quadrant((0.0, 3.6),   (4.5,  100), 1),
            TSO.Quadrant((3.6, 4.5),   (-100, 1.0), 4),
            TSO.Quadrant((4.5, 100.0), (-100, 1.0), 3),
            TSO.Quadrant((3.6, 100.0), (1.0,  100), 1)
		], =#=#
        maxiter=5000, display=:none)


# ╔═╡ 7ea0d77a-d854-490d-94fc-da997a00649a
fromT_toE(new_eos_folder, new_eosE_folder, upsample=2048)

# ╔═╡ Cell order:
# ╟─317e3358-8512-4166-9596-23ae7bb380c6
# ╠═3c3d9886-587b-4209-a421-1f21e6749e04
# ╟─48efb4ac-ca3f-4a48-b5ed-9e173fb01846
# ╟─2de77db3-f110-4e20-b68e-b00ce9451e46
# ╠═aad8b35a-ea58-11ed-2350-79fec9c7f69c
# ╟─9aae590d-b656-43f3-949e-0966f1a7a32a
# ╟─9fb3f203-331e-4b76-b94c-c992f625c80a
# ╠═b45ecd12-4a70-42ac-bbee-9d6568adc214
# ╟─0162409f-700e-457f-bd2d-ca685a025274
# ╠═94c6f318-22be-47e9-9248-adf5877e4d8a
# ╠═a8c93f35-2ee0-4e62-b9f2-4762d108f02e
# ╠═7ea0d77a-d854-490d-94fc-da997a00649a
