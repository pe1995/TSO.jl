"""
    compute_formation_opacities(table_folder, av_path, name=""; logg=log10(2.75e4))

Compute the rosseland optical depth for the entire table. Afterwards, compute the formation height of every wavelength point in the unbinned opacity table.
Note: `name` will decide the name suffix under which the formation opacities will be saved. The path will be (within the `table_folder`)
`combined_formation_opacities_name_extension.hdf5`. The `extension` is assumed to be a short-hand for the given unbinned table.
If under this `extension` at the very end no table with the name `combined_formation_opacities_name_extension.hdf5` is found, it will be created.
Also a new opacity table with updated rosseland opacities is saved.
"""
function compute_formation_opacities(name, table_folder, eos_path, opa_path, av_path; corr_χ=nothing, corr_S=nothing, logg=log10(2.75e4))
    opacities = reload(
        SqOpacity, 
		joinpath(table_folder, opa_path), 
        mmap=false
    )
    
    eos = reload(
        SqEoS,     
		joinpath(table_folder, eos_path)
    )

    if !isnothing(corr_χ)
        @info "Correcting formation opacities."
        TSO.scale!(opacities; κ=corr_χ, src=corr_S)
    end
    
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
                                        corr_χ=nothing, corr_S=nothing, 
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
	
    weights = ω_midpoint(opacities)
    model = @optical flip(Average3D(eos, av_path, logg=logg), depth=true) eos opacities
    
    bin = if Nbins > 1
        formOpacities = reload(SqOpacity, 
            joinpath(table_folder, formation_opacities), 
            mmap=true
        )

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
    else
        # in the monochromatic case we dont need formation opacities and dont need binning
        ones(length(opacities.λ))
    end

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
	
    binned_opacities = tabulate(bin, weights, eos, opacities, sopacities, transition_model=model; corr_χ=corr_χ, corr_S=corr_S)

    binned_opacities = if Nbins == 1
        @info "Replacing binned opacity with Rosseland opacity."

        # Replace binned opacity with rosseland opacity
        binned_opacities.opacities.κ[:, :, 1] .= binned_opacities.opacities.κ_ross[:, :]

        # multiply with density from EoS
        TSO.@binned binned_opacities.opacities eos 
    else
        binned_opacities
    end
    
    save_table(binned_opacities)
end



"""
    convert_fromT_toE(table_folder, folder_new)

Convert the binned opacities + eos from a T-ρ to a Eint-ρ grid. Upsample the resulting table.
"""
function convert_fromT_toE(table_folder, folder_new; upsample=1000, extend=false, downD=0.4, downE=0.1, upD=0.01, upE=0.01, lnEimin=nothing, lnEimax_V=nothing, lnEimin_V=nothing, lnEimax=nothing, eos_radius=-1, opa_radius=-1)
    eos = reload(SqEoS, joinpath(table_folder, "eos.hdf5"))
    opa = reload(SqOpacity, joinpath(table_folder, "binned_opacities.hdf5"))
    aos = @axed eos

    # make sure new dir exists
    !isdir(folder_new) && mkdir(folder_new) 

    # Before switching energy sources, make sure that everything is monotonic
    TSO.smoothAccumulate!(aos, spline=true)

    eos, opa = if (eos_radius>0) & (opa_radius>0)
        @info "Smoothing EoS + opacities with Gaussian kernel radii $(eos_radius) (EoS) and $(opa_radius) (opacities)"
        TSO.smooth(eos, opa, eos_radius=eos_radius, opa_radius=opa_radius)
    else
        eos, opa
    end

     # save the T-grid EoS to make comparison easier later
     save(opa, joinpath(folder_new, "binned_opacities_T.hdf5"))
     save(eos, joinpath(folder_new, "eos_T.hdf5"))

     # also copy over bin assignment, if present
     if isfile(joinpath(table_folder, "bin_assignment.hdf5"))
        cp(joinpath(table_folder, "bin_assignment.hdf5"), joinpath(folder_new, "bin_assignment.hdf5"), force=true)
    end

    # also save per unit volume
    eos_ev = deepcopy(eos)
	for j in axes(eos_ev.lnEi, 2)
		eos_ev.lnEi[:, j] .+= eos_ev.lnRho[j]
	end
    aos_ev = @axed(eos_ev)

    # regular table
    _convert_fromT_toE(
        aos, opa, folder_new, 
        upsample, extend, 
        downD, downE, upD, upE, 
        lnEimin, lnEimax, eos_radius, 
        opa_radius; name=""
    )

    # ev table
    _convert_fromT_toE(
        aos_ev, opa, folder_new, 
        upsample, extend, 
        downD, downE, upD, upE, 
        lnEimin_V, lnEimax_V, eos_radius, 
        opa_radius; name="_V"
    )
end

_convert_fromT_toE(aos, opa, folder_new,
    upsample=1000, extend=false, 
    downD=0.4, downE=0.1, upD=0.01, upE=0.01, 
    lnEimin=nothing, lnEimax=nothing, 
    eos_radius=-1, opa_radius=-1; name="") = begin

    # possibly extrapolate EoS
    eos_new, opa_new = if extend
        @info "Extrapolating EoS beyond limits."
        TSO.extend(aos, opa, downD=downD, downE=downE, upD=upD, upE=upE)
    else
        aos.eos, opa
    end

    # Again making it monotonic
    TSO.smoothAccumulate!(@axed(eos_new), spline=true)

    # Limit the internal energy range to what is occupied by the initial model to make sure
    # the model atmosphere will be well sampled
    if !isnothing(lnEimax)
        lnEimax_orig = maximum(eos_new.lnEi)
        mask = eos_new.lnEi.>=lnEimax
        eos_new.lnEi[mask] .= lnEimax

        if !isnothing(lnEimin)
            @info "Limiting the internal energy (log) between $(lnEimin) - $(lnEimax) on the energy grid." 
            @info "The original limits were $(minimum(eos_new.lnEi)) - $(lnEimax_orig)." 
            eos_new.lnEi[eos_new.lnEi.<=lnEimin] .= lnEimin
        else
            @info "Limiting the internal energy (log) to $(lnEimax) on the energy grid." 
            @info "The original limits was $(lnEimax_orig)." 
        end
        @info "The upper limit will have an effect on $(count(mask)/length(mask)*100)% of points in the rho-T table."
    end

    #@show log10(exp(minimum(eos_new.lnEi))) log10(exp(maximum(eos_new.lnEi)))
    eosE, opaE = switch_energy(@axed(eos_new), opa_new, upsample=upsample, conservative=false);
    aosE = @axed eosE
    #@show log10(exp(minimum(eosE.lnEi))) log10(exp(maximum(eosE.lnEi)))

    # do not allow values larger than in the original table
    for (i, l) in enumerate(opaE.λ)
        limit_beyond!(view(opaE.κ, :, :, i), view(opa_new.κ, :, :, i))
        limit_beyond!(view(opaE.src, :, :, i), view(opa_new.src, :, :, i))
    end
    limit_beyond!(eosE.lnPg, eos_new.lnPg)
    limit_beyond!(eosE.lnNe, eos_new.lnNe)
    limit_beyond!(eosE.lnRoss, eos_new.lnRoss)
    limit_beyond!(eosE.lnT, eos_new.lnT)
    
    TSO.fill_nan!(aosE, opaE)
    #TSO.set_limits!(aosE, opaE)

    # artificially remove sharp edges in the table at low density (test)
    #=tlim = lookup(eosE, :lnEi, eosE.lnRho, fill!(similar(eosE.lnRho), log(4500)))
	κ_limit = [lookup(eosE, opaE, :κ, eosE.lnRho, tlim, b) for b in eachindex(opaE.λ)]
	rho_limit = 1e-1
	
	κ_limited = deepcopy(opaE.κ)	
	for j in axes(κ_limited, 2)
		for i in axes(κ_limited, 1)
			if (eosE.lnEi[i] <= tlim[j]) & ((eosE.lnRho[j] |> exp) < rho_limit)
				κ_limited[i, j, :] .= [κ_limit[b][j] for b in eachindex(opaE.λ)]
			end
		end
	end
    opaE.κ .= κ_limited=#

    for_dispatch(eosE, opaE.κ, opaE.src, ones(eltype(opaE.src), size(opaE.src)...), folder_new, name=name)

    save(opaE, joinpath(folder_new, "binned_opacities$(name).hdf5"))
    save(eosE, joinpath(folder_new, "eos$(name).hdf5"))
end


function limit_beyond!(new_arr, old_arr; set_to=NaN)
    m = (new_arr .> maximum(old_arr)) .|| (new_arr .< minimum(old_arr))
    new_arr[m] .= set_to
end

"""
    convert_fromT_toE(table_folder, folder_new, av_path; lnEi_stretch=1.0, kwargs...)

Convert the binned opacities + eos from a T-ρ to a Eint-ρ grid. If needed, 
the internal energy is cut where the models internal energy ends (+lnEi_stretch*absolute difference top-bottom).
"""
convert_fromT_toE(table_folder, folder_new, av_path; lnEi_stretch=1.0, kwargs...) = begin
    eos = reload(SqEoS, joinpath(table_folder, "eos.hdf5"))
    opa = reload(SqOpacity, joinpath(table_folder, "binned_opacities.hdf5"));
    lnEimin, lnEilimit = _get_e_limit(eos, opa, av_path, lnEi_stretch)

    # also save per unit volume
    eos_ev = deepcopy(eos)
	for j in axes(eos_ev.lnEi, 2)
		eos_ev.lnEi[:, j] .+= eos_ev.lnRho[j]
	end
    lnEimin_V, lnEilimit_V = _get_e_limit(eos_ev, opa, av_path, lnEi_stretch)

    convert_fromT_toE(table_folder, folder_new; lnEimin=lnEimin, lnEimax=lnEilimit, lnEimin_V=lnEimin_V, lnEimax_V=lnEilimit_V, kwargs...)
end

_get_e_limit(eos, opa, av_path, lnEi_stretch) = begin
    model = @optical flip(Average3D(eos, av_path), depth=true) eos opa

    lnEimin = minimum(model.lnEi)
    lnEimax = maximum(model.lnEi)

    # determin the temperature minimum
    tmin = minimum(eos.lnT)

    # lookup the corresponding internal energy curve
    lnEi_min_curve = eos.lnEi[1, :]

    # use the median to avoid overly small lower limits
    lnEimin = median(lnEi_min_curve)

    (lnEimin, lnEimax + abs(lnEimax - lnEimin) * lnEi_stretch)
end


function create_E_from_T(table_folder, name="", av_path=nothing; 
                                upsample=2048,
                                name_extension="DIS_MARCS_E",
                                version="v0.1", lnEi_stretch=1.0)
    new_table_name = TSO.join_full(name_extension, name, version, add_start=false)

    if !isnothing(av_path)
        convert_fromT_toE(table_folder, new_table_name, av_path, upsample=upsample)
    else
        convert_fromT_toE(table_folder, new_table_name, upsample=upsample)
    end
end