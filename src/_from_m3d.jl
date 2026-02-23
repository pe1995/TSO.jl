# ============================================================================
# Input for M3D opacity table computation
# ============================================================================

function eosTableInput(path; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3)
	if isdir(path)
		# clear the folder
		rm(path, recursive=true)
	end
	mkdir(path)

	z = [-1.0, 1.0]
	saveAsText(
		joinpath(path, "TSO-M3D"), 
		z=z, 
		T=[minT, maxT], 
		ρ=[minρ, maxρ]
	)

	"TSO-M3D"
end

function saveAsText(path_new; z, T, ρ, vmic=1.0)
	k = length(z)
	p = zeros(eltype(T), k)

	if first(ρ) > last(ρ)
		reverse!(z)
		reverse!(T)
		reverse!(ρ)
	end
	if first(z) > last(z)
		z .*= -1
	end
	
	open(path_new, "w") do f
		write(f, "TSO-M3D-table-column\n")
		write(f, "$(k)\n")
		for i in eachindex(z)
			write(f, "$(z[i]) $(T[i]) $(p[i]) $(ρ[i]) $(vmic)")
			if i < length(z)
				write(f, "\n")
			end
		end
	end

	path_new
end

function EoSTableInput(path; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3, vmic=1.0, outputname="TSO-M3D")
	if isdir(path)
		rm(path, recursive=true)
	end
	mkdir(path)

	z = [-1.0, 1.0]
	TSO.saveAsText(
		joinpath(path, outputname), 
		z=z, 
		T=[minT, maxT], 
		ρ=[minρ, maxρ],
		vmic=vmic
	)

	outputname
end

"""
    computeEosTable(model; folder, linelist, λ_file, λs, λe, δλ, δlnT, δlnρ, FeH, nν,
				in_log=true, slurm=false, m3dis_kwargs=Dict(), abund_file, tmolim, kwargs...)

Compute opacity and EoS tables using M3D.
"""
function computeEoSTable(model; folder, linelist, λ_file, λs, λe, δλ, δlnT, δlnρ, FeH, nν,
				in_log=true, slurm=false, m3dis_kwargs=Dict(), abund_file, tmolim, kwargs...)
	spec_para = (!isnothing(λ_file)) ? Dict(:lam_file=>λ_file, :in_air=>false) : Dict(:daa=>δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>in_log, :in_air=>false)
    MUST.whole_spectrum(
		model, 
		namelist_kwargs=(
			:model_folder=>folder,
			:linelist=>nothing,
			:absmet=>nothing,
			:linelist_params=>(:linelist_folder=>linelist,),
			:atom_params=>(:atom_file=>"", ),
			:spectrum_params=>(spec_para...,),
			:atmos_params=>(
				:dims=>1, 
				:atmos_format=>"Text",
				:use_rho=>true, 
				:use_ne=>false,
				:FeH=>FeH,
				:nz=>2,
				:amr=>false
			),
			:m3d_params=>(
				:n_nu=>nν, 
				:ilambd=>0,
				:short_scheme=>"disk_center",
				:long_scheme=>"none",
				:make_eos=>true
			),
			:composition_params=>(
                :abund_file=>abund_file,
				:ldtemp=>δlnT,
				:ldrho=>δlnρ,
				:tmolim=>tmolim,
				:mhd_eos=>false,
			),
            kwargs...
		),
		m3dis_kwargs=m3dis_kwargs,
		slurm=slurm
	)
end

# ============================================================================
# Get Opacity and EoS tables from TUMULT
# ============================================================================

function _get_opacity(run, from)
	chi = pyconvert(
		Array, 
		run.run.read_patch_save(from, concat=true, fdim=0, lazy=false)[0]
	)

	chi  = chi[:, :, :] 
	l = Base.convert.(eltype(chi), reverse(pyconvert(Array, numpy.fromfile(joinpath("$(run.run.sfolder)", "out_lam.bin"), dtype=numpy.float64))))
	if !issorted(l)
		@warn "Wavelength array appears to be not sorted. This may indicate a bug."
		m = sortperm(l)
		l .= l[m]
		chi .= chi[m, :, :]
	end

	chi[chi .< 1e-30] .= NaN
	chi[chi .> 1e30] .=  NaN

	return l, permutedims(chi, (2, 3, 1))
end

function _get_opacity_lazy(run, from)
	chi_lazy = run.run.read_patch_save(from, concat=true, fdim=0, lazy=true)[0]
	
	shape = pyconvert(Tuple, chi_lazy.shape)
	if pytruth(pybuiltins.hasattr(chi_lazy, "dtype"))
		dtype_str = pyconvert(String, pybuiltins.str(chi_lazy.dtype))
	else
		dtype_str = pyconvert(String, pybuiltins.str(chi_lazy.mlist[0].dtype))
	end
	T = dtype_str == "float64" ? Float64 : Float32

	tmpfile = tempname() * ".bin"
	
	# The original array has shape (N_lam, N_T, N_rho)
	# The output needs to be permuted to (N_T, N_rho, N_lam)
	permuted_shape = (shape[2], shape[3], shape[1])
	chi_mmap = Mmap.mmap(open(tmpfile, "w+"), Array{T, length(permuted_shape)}, permuted_shape)

	if pytruth(pybuiltins.hasattr(chi_lazy, "mlist"))
		caxis = pyconvert(Int, chi_lazy.caxis) + 1 # Julia is 1-indexed
		cstarts = pyconvert(Vector{Int}, chi_lazy.cstart)
		
		mlist = pybuiltins.list(chi_lazy.mlist)
		for (i, m) in enumerate(mlist)
			m_array = pyconvert(Array, m)
			
			start_idx = cstarts[i] + 1
			end_idx = start_idx + size(m_array, caxis) - 1
			
			indices_raw = Any[Colon() for _ in 1:length(shape)]
			indices_raw[caxis] = start_idx:end_idx
			
			permuted_chunk = permutedims(m_array, (2, 3, 1))
			permuted_caxis = caxis == 1 ? 3 : (caxis == 2 ? 1 : 2)
			indices = Any[Colon() for _ in 1:length(permuted_shape)]
			indices[permuted_caxis] = start_idx:end_idx
			
			chi_mmap[indices...] .= permuted_chunk
		end
	else
		# Plain array/memmap, copy in chunks along the first dimension (lam)
		chunk_size = max(1, 100_000_000 ÷ (prod(shape[2:end])))
		for i in 1:chunk_size:shape[1]
			end_idx = min(i + chunk_size - 1, shape[1])
			
			chunk = pyconvert(Array, chi_lazy[i-1:end_idx-1]) # python is 0-indexed
			permuted_chunk = permutedims(chunk, (2, 3, 1))
			
			chi_mmap[.., i:end_idx] .= permuted_chunk
		end
	end

	l = Base.convert.(T, reverse(pyconvert(Array, numpy.fromfile(joinpath("$(run.run.sfolder)", "out_lam.bin"), dtype=numpy.float64))))
	
	if !issorted(l)
		@warn "Wavelength array appears to be not sorted. This may indicate a bug."
		m = sortperm(l)
		l .= l[m]
		chi_mmap .= chi_mmap[:, :, m]
	end

	for i in eachindex(chi_mmap)
		val = chi_mmap[i]
		if val < 1e-30 || val > 1e30
			chi_mmap[i] = NaN
		end
	end

	return l, chi_mmap
end

"""
    get_opacity(run, from="opa"; mmap=false)

Extract EoS + Opacities from given run. Return arrays.
If `mmap=true`, it writes Python chunks to a Julia-mapped file 
instead of holding everything in RAM, which is more memory-efficient.
"""
function get_opacity(run, from="opa"; mmap=false)
	if !numpy_loaded[]
        load_numpy!()
    end

    if !mmap
        return _get_opacity(run, from)
    else
        return _get_opacity_lazy(run, from)
    end
end

"""
    get_opacity(run)

Extract EoS + Opacities from given run. Return arrays
"""
function get_eos(run)
    temp = pyconvert(Array, run.run.temp)[:, 1]
    rho  = pyconvert(Array, run.run.rho)[:, :]
    ne   = pyconvert(Array, run.run.ne)[:, :]
    pg   = pyconvert(Array, run.run.pg)[:, :]
    E    = pyconvert(Array, run.run.E)[:, :] .+ Base.convert(eltype(rho), (13.595+5.0)/2.380491e-24*1.60218e-12)
    χ500 = pyconvert(Array, run.run.chi)[:, :]
	
	#=eos = pyconvert(
		Array, 
		run.run.read_patch_save("eos", concat=false)[1,1]
	)

    temp = eos[1, :, 1]
    rho  = eos[2, 1, :]
    ne   = eos[3, :, :]
    pg   = eos[4, :, :]
    E    = eos[5, :, :] .+ Base.convert(eltype(rho), (13.595+5.0)/2.380491e-24*1.60218e-12)=#
	
	E[E .< 1e-30] .= NaN
	E[E .> 1e30] .= NaN

	pg[pg .< 1e-30] .= NaN
	pg[pg .> 1e30] .= NaN

	ne[ne .< 1e-30] .= NaN
	ne[ne .> 1e30] .= NaN

	χ500[χ500 .< 1e-30] .= NaN
	χ500[χ500 .> 1e30] .= NaN
	
	log.(rho[1, :]), log.(temp), log.(E), log.(pg), log.(ne), log.(χ500 ./ rho)
end

# ============================================================================
# Add Opacities from multiple tables
# ============================================================================

function add_opacities!(opa1, opa...)
	for i in eachindex(opa)
		opa1.κ .+= opa[i].κ
	end
end

# ============================================================================
# Convert Opacity and EoS from TUMULT to RegularOpacityTable
# ============================================================================

"""
    collect_opacity(run; compute_ross=false, mini=false, mmap=false)

Extract EoS + Opacities from given run. Collect them together in a table.
If `mini=true`, it skips memory allocation for the Source function array `S`.
If `mmap=true`, it uses file memory-mapping for parsing the NumPy arrays and creating `S`.
"""
function collect_opacity(run; compute_ross=true, mini=false, mmap=true)
	rho, temp, E, pg, ne, χ500 = get_eos(run)
	l, chi = get_opacity(run, "opa"; mmap=mmap)
	_, scat = try
		get_opacity(run, "scat"; mmap=mmap)
	catch
		@warn "No scattering opacity found."
		nothing, nothing
	end

	ross = fill!(similar(pg), 1.0)
    
    S = if mini
        fill!(similar(pg, 1, 1, 1), 1.0)
    elseif mmap
        tmpfile_S = tempname() * ".bin"
        S = Mmap.mmap(open(tmpfile_S, "w+"), Array{eltype(chi), ndims(chi)}, size(chi))
        fill!(S, 1.0)
		S
    else
	    fill!(similar(chi), 1.0)
    end

	eos = TSO.SqEoS(
		rho, temp, E, pg, log.(ross), ne
	)
	eos500 = TSO.SqEoS(
		rho, temp, E, pg, χ500, ne
	)
    aos = @axed eos
    aos500 = @axed eos500


	chi = TSO.SqOpacity(
		chi, ross, S, l, false
	)
	scat = if !isnothing(scat)
		TSO.SqOpacity(
			scat, ross, S, l, false
		)
	else
		nothing
	end

	# add scattering opacity to normal opacity
	#if !isnothing(scat)
	#	@info "Adding scattering opacity to absorption."
	#	add_opacities!(chi, scat)
	#end

	# compute rosseland opacity
	add_radiation_quantities!(eos, chi, scat, compute_ross=compute_ross)
	fill_nan!(aos, chi, scat)
	fill_nan!(aos500)
	set_limits!(aos, chi, small=1e-30, large=1e30)
	set_limits!(aos, scat, small=1e-30, large=1e30)
	set_limits!(aos500, small=1e-30, large=1e30)

	eos, eos500, chi, scat
end

function add_radiation_quantities!(eos, opa, scat=nothing; compute_ross=false)
	if size(opa.src) != (1, 1, 1)
		T = exp.(eos.lnT)
		for k in axes(opa.src, 3)
			for j in axes(opa.src, 2)
				for i in axes(opa.src, 1)
					opa.src[i, j, k] = TSO.Bλ_fast(opa.λ[k], T[i])
				end
			end
		end
	end

	if compute_ross
		@info "Computing table Rosseland opacity. This may take a few minutes."
		rosseland_opacity!(eos.lnRoss, @axed(eos), opa)
		transfer_rosseland!(@axed(eos), opa)
		if !isnothing(scat)
			transfer_rosseland!(@axed(eos), scat)
		end
	end

	eos, opa, scat
end