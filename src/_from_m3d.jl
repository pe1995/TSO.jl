#= Running M3D =#

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










#= Converting the output =#

"""
    get_opacity(run)

Extract EoS + Opacities from given run. Return arrays
"""
function get_opacity(run, from="opa")
	chi = pyconvert(
		Array, 
		run.run.read_patch_save(from, concat=true, fdim=0, lazy=false)[0]
	)

	chi  = chi[:, :, :] 
	l = Base.convert.(eltype(chi), run.lam)

	chi[chi .< 1e-30] .= NaN
	chi[chi .> 1e30] .=  NaN

	l, permutedims(chi, (2, 3, 1))
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

function add_opacities!(opa1, opa...)
	for i in eachindex(opa)
		opa1.κ .+= opa[i].κ
	end
end




"""
    collect_opacity(run)

Extract EoS + Opacities from given run. Collect them together in a table.
"""
function collect_opacity(run; compute_ross=false)
	rho, temp, E, pg, ne, χ500 = get_eos(run)
	l, chi = get_opacity(run, "opa")
	_, scat = try
		get_opacity(run, "scat")
	catch
		@warn "No scattering opacity found."
		nothing, nothing
	end

	ross = fill!(similar(pg), 1.0)
	S = fill!(similar(chi), 1.0)

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
	if !isnothing(scat)
		@info "Adding scattering opacity to absorption."
		add_opacities!(chi, scat)
	end

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
	T = exp.(eos.lnT)
	for j in eachindex(eos.lnRho)
		opa.src[:, j, :] .= TSO.Bλ(opa.λ, T)
	end

	if compute_ross
		@info "Computing table Rosseland opacity. This may take a while..."
		rosseland_opacity!(eos.lnRoss, @axed(eos), opa)
		transfer_rosseland!(@axed(eos), opa)
		if !isnothing(scat)
			transfer_rosseland!(@axed(eos), scat)
		end
	end

	eos, opa, scat
end