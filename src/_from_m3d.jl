# ============================================================================
# Input for M3D opacity table computation
# ============================================================================

function eosTableInput(path; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3)
	if isdir(path)
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

# ============================================================================
# Get Opacity and EoS tables from TUMULT
# ============================================================================

function _get_opacity(run, from)
	chi = pyconvert(
		Array,
		run.run.read_patch_save(from, concat=true, fdim=0, lazy=false)[0]
	)

	chi  = chi[:, :, :]
	l = Base.convert.(eltype(chi), pyconvert(Array, numpy.fromfile(joinpath("$(run.run.sfolder)", "out_lam.bin"), dtype=numpy.float64)))
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

	tmpfile = tempname(pwd(), suffix=".bin")

	permuted_shape = (shape[2], shape[3], shape[1])
	chi_mmap = Mmap.mmap(open(tmpfile, "w+"), Array{T, length(permuted_shape)}, permuted_shape)

	if pytruth(pybuiltins.hasattr(chi_lazy, "mlist"))
		caxis = pyconvert(Int, chi_lazy.caxis) + 1
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
		chunk_size = max(1, 100_000_000 ÷ (prod(shape[2:end])))
		for i in 1:chunk_size:shape[1]
			end_idx = min(i + chunk_size - 1, shape[1])

			chunk = pyconvert(Array, chi_lazy[i-1:end_idx-1])
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

function get_eos(run)
    temp = pyconvert(Array, run.run.temp)[:, 1]
    rho  = pyconvert(Array, run.run.rho)[:, :]
    ne   = pyconvert(Array, run.run.ne)[:, :]
    pg   = pyconvert(Array, run.run.pg)[:, :]
    E    = pyconvert(Array, run.run.E)[:, :] .+ Base.convert(eltype(rho), (13.595+5.0)/2.380491e-24*1.60218e-12)
    χ500 = pyconvert(Array, run.run.chi)[:, :]

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
# Convert Opacity and EoS from TUMULT to EquationOfState
# ============================================================================

function collect_opacity(run; compute_ross=true, mini=false, mmap=true)
	rho, temp, E, pg, ne, χ500 = get_eos(run)
	l, chi = get_opacity(run, "opa"; mmap=mmap)
	_, scat = try
		get_opacity(run, "scat"; mmap=mmap)
	catch
		@warn "No scattering opacity found."
		nothing, nothing
	end

	assemble_opacity(rho, temp, E, pg, ne, χ500, l, chi, scat; compute_ross=compute_ross, mini=mini, mmap=mmap)
end

"""
    assemble_opacity(rho, temp, E, pg, ne, χ500, l, chi, scat; compute_ross=true, mini=false, mmap=true)

Monochromatic table from the M3D output. NaNs are filled (`nan_mask`, `nan_mask500` mark
where), the source function is the Planck function and the Rosseland opacity is integrated.
With `mini=true` no source function is stored.
"""
function assemble_opacity(rho, temp, E, pg, ne, χ500, l, chi, scat; compute_ross=true, mini=false, mmap=true)
	ross = fill!(similar(pg), 1.0)

    S = if mini
        nothing
    elseif mmap
        tmpfile_S = tempname(pwd(), suffix=".bin")
        S = Mmap.mmap(open(tmpfile_S, "w+"), Array{eltype(chi), ndims(chi)}, size(chi))
        fill!(S, 1.0)
		S
    else
	    fill!(similar(chi), 1.0)
    end

	tab = EquationOfState(rho, temp, E, pg, log.(ross), ne;
		opacity=OpacityData(chi, S, l, false),
		scattering=isnothing(scat) ? nothing : OpacityData(scat, S, l, false),
		meta=Dict{Symbol, Any}(:kind => "monochromatic")
	)
	tab.eos[:lnK500] = χ500

	tab.eos[:nan_mask] = fill_nan!(tab; eos_fields=(:lnEi, :lnPg, :lnNe, :lnRoss))
	tab.eos[:nan_mask500] = fill_nan!(tab; eos_fields=(:lnEi, :lnPg, :lnNe, :lnK500), opacities=false)

	set_limits!(tab, small=1e-30, large=1e30)

	add_radiation_quantities!(tab, compute_ross=compute_ross)

	tab
end

function add_radiation_quantities!(tab::EquationOfState; compute_ross=false)
	if haskey(tab.opacity, :src)
		_planck_source!(tab.opacity[:src], tab.opacity.λ, exp.(tab.lnT))
	end

	if compute_ross
		rosseland_opacity!(tab.eos[:lnRoss], tab)
	end

	tab
end

function _planck_source!(src, λ, T)
	for k in axes(src, 3)
		for j in axes(src, 2)
			for i in axes(src, 1)
				src[i, j, k] = TSO.Bλ_fast(λ[k], T[i])
			end
		end
	end
end
