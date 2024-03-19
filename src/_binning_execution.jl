"""
    compute_formation_opacities(table_folder, av_path, name=""; logg=log10(2.75e4))

Compute the rosseland optical depth for the entire table. Afterwards, compute the formation height of every wavelength point in the unbinned opacity table.
Note: `name` will decide the name suffix under which the formation opacities will be saved. The path will be (within the `table_folder`)
`combined_formation_opacities_name_extension.hdf5`. The `extension` is assumed to be a short-hand for the given unbinned table.
If under this `extension` at the very end no table with the name `combined_formation_opacities_name_extension.hdf5` is found, it will be created.
Also a new opacity table with updated rosseland opacities is saved.
"""
function compute_formation_opacities(name, table_folder, eos_path, opa_path, av_path; logg=log10(2.75e4))
    opacities = reload(
        SqOpacity, 
		joinpath(table_folder, opa_path), 
        mmap=true
    )
    
    eos = reload(
        SqEoS,     
		joinpath(table_folder, eos_path)
    )
    
    aos = @axed eos
    model = flip(Average3D(av_path, logg=logg), depth=true)
    τ_ross, τ_λ = optical_depth(aos, opacities, model)
    d_ross, d_κ = formation_height(model, aos, opacities, τ_ross, τ_λ)
    formation_opacities = SqOpacity(d_κ, d_ross, opacities.src, opacities.λ, true)
	
    save(
        formation_opacities, 
		joinpath(table_folder, "combined_formation_opacities_$(name).hdf5")
    )
end


function compute_rosseland_opacities(table_folder, eos_path, opa_path)
    opacities = reload(
        SqOpacity, 
		joinpath(table_folder, opa_path), 
        mmap=false
    )
    
    eos = reload(
        SqEoS,     
		joinpath(table_folder, eos_path)
    )
    
    aos = @axed eos
  
    @info "computing rosseland"
    rosseland_opacity!(aos.eos.lnRoss, aos, opacities; 
                weights=ω_midpoint(opacities))
		
    transfer_rosseland!(aos.eos, opacities)
		
    eos_ross_name = "ross_"*eos_path
    save(aos.eos,   joinpath(table_folder, eos_ross_name))
    save(opacities, joinpath(table_folder, opa_path))

    eos_ross_name
end



"""
    bin_opacity_table(table_folder, av_path, new_folder; logg=log10(2.75e4))

Bin opacities, based on a specific binning method, using formation opacities computed previously. Save the binned tables.
"""
function bin_opacity_table(name, table_folder, eos_path, opa_path, av_path; 
                                        method=:uniform, 
                                        use_contribution=false, 
                                        logg=log10(2.75e4), 
                                        Nbins=5, 
                                        line_bins=3, 
                                        scattering_path=nothing,
                                        name_extension="DIS_MARCS",
                                        formation_opacities="combined_formation_opacities_$(name).hdf5",
                                        version="v0.1",
                                        kwargs...)
    eos = reload(
        SqEoS, 
		joinpath(table_folder, eos_path)
    ) 
	
    opacities = reload(
        SqOpacity, 
		joinpath(table_folder, opa_path), 
		mmap=true
    )

    sopacities = isnothing(scattering_path) ? nothing : reload(
        SqOpacity, 
		joinpath(table_folder, scattering_path), 
		mmap=true
    )
	
    formOpacities = reload(SqOpacity, 
		joinpath(table_folder, formation_opacities), 
		mmap=true
    )

    weights = ω_midpoint(opacities)
    model = @optical flip(Average3D(eos, av_path, logg=logg), depth=true) eos opacities

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

    # saving the table
    save_table(binned_opacities) = begin
        eos_table_name = TSO.join_full(name_extension, name, version, add_start=false)

        !isdir(eos_table_name) && mkdir(eos_table_name) 

        cp(joinpath(table_folder, eos_path), 
			joinpath(eos_table_name, "eos.hdf5"), force=true)

        save(binned_opacities.opacities, joinpath(eos_table_name, "binned_opacities.hdf5"))

        fid = TSO.HDF5.h5open(joinpath(eos_table_name, "bin_assignment.hdf5"), "w")
        fid["bins"] = bin
        fid["lambda"] = opacities.λ
        close(fid)

        eos_table_name
    end
	
    binned_opacities = if isnothing(sopacities)
        tabulate(bin, weights, eos, opacities, transition_model=model)
    else
        tabulate(bin, weights, eos, opacities, sopacities, transition_model=model)
    end
    
    save_table(binned_opacities)
end



"""
    convert_fromT_toE(table_folder, folder_new)

Convert the binned opacities + eos from a T-ρ to a Eint-ρ grid. Upsample the resulting table.
"""
function convert_fromT_toE(table_folder, folder_new; upsample=1000, extend=false, downD=0.4, downE=0.1, upD=0.01, upE=0.01)
    eos = reload(SqEoS,     joinpath(table_folder, "eos.hdf5"))
    opa = reload(SqOpacity, joinpath(table_folder, "binned_opacities.hdf5"));
    aos = @axed eos

    # Before switching energy sources, make sure that everything is monotonic
    TSO.smoothAccumulate!(aos, spline=true)
    eos_new, opa_new = if extend
        @info "Extrapolating EoS at $(table_folder) beyond limits."
        TSO.extend(aos, opa, downD=downD, downE=downE, upD=upD, upE=upE)
    else
        aos.eos, opa
    end
    TSO.smoothAccumulate!(@axed(eos_new), spline=true)

    eosE, opaE = switch_energy(@axed(eos_new), opa_new, upsample=upsample, conservative=false);
    aosE = @axed eosE

    TSO.fill_nan!(aosE, opaE)
    TSO.set_limits!(aosE, opaE)

    !isdir(folder_new) && mkdir(folder_new) 

    for_dispatch(eosE, opaE.κ, opaE.src, ones(eltype(opaE.src), size(opaE.src)...), folder_new)

    save(opaE, joinpath(folder_new, "binned_opacities.hdf5"))
    save(eosE, joinpath(folder_new, "eos.hdf5"))

    # also save the T-grid EoS to make comparison easier later
    save(opa, joinpath(folder_new, "binned_opacities_T.hdf5"))
    save(eos, joinpath(folder_new, "eos_T.hdf5"))
end

function create_E_from_T(table_folder, name=""; 
                                upsample=2048,
                                name_extension="DIS_MARCS_E",
                                version="v0.1")
    new_table_name = TSO.join_full(name_extension, name, version, add_start=false)
    convert_fromT_toE(table_folder, new_table_name, upsample=upsample)
end






bin_assignment(path) = begin
	o = TSO.HDF5.h5open(path)
	TSO.HDF5.read(o["bins"])
end