#= Running M3D =#

function eosTableInput(path; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3)
	if isdir(path)
		# clear the folder
		rm(path, recursive=true)
	end
	mkdir(path)

	z = [-1.0, 0.0, 1.0]
	saveAsText(
		joinpath(path, "TSO-M3D"), 
		z=z, 
		T=[minT, (maxT-minT)/2, maxT], 
		ρ=[minρ, (maxρ-minρ)/2, maxρ]
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
function get_opacity(run)
	chi = pyconvert(
		Array, 
		run.run.read_patch_save("opa", concat=true, fdim=0, lazy=false)[0]
	)
	#=eos = pyconvert(
		Array, 
		run.run.read_patch_save("eos", concat=true, fdim=3, lazy=false)[0]
	)=#
	eos = pyconvert(
		Array, 
		run.run.read_patch_save("eos", concat=false)[1,1]
	)

	chi  = chi[:, :, :]
    temp = eos[1, :, 1]
    rho  = eos[2, 1, :]
    ne   = eos[3, :, :]
    pg   = eos[4, :, :]
    E    = eos[5, :, :] .+ Base.convert(eltype(rho), (13.595+5.0)/2.380491e-24*1.60218e-12)
	
	l = Base.convert.(eltype(rho), run.lam)

	chi[chi .< 1e-30] .= NaN
	chi[chi .> 1e30] .=  NaN
	
	E[E .< 1e-30] .= NaN
	E[E .> 1e30] .= NaN

	pg[pg .< 1e-30] .= NaN
	pg[pg .> 1e30] .= NaN

	ne[ne .< 1e-30] .= NaN
	ne[ne .> 1e30] .= NaN
	
	l, log.(rho), log.(temp), log.(E), log.(pg), log.(ne), permutedims(chi, (2, 3, 1))
end

"""
    collect_opacity(run)

Extract EoS + Opacities from given run. Collect them together in a table.
"""
function collect_opacity(run; compute_ross=false)
	l, rho, temp, E, pg, ne, chi = get_opacity(run)

	ross = fill!(similar(pg), 1.0)
	S = fill!(similar(chi), 1.0)

	eos = TSO.SqEoS(
		rho, temp, E, pg, log.(ross), ne
	)

	opa = TSO.SqOpacity(
		chi, ross, S, l, false
	)

    aos = @axed eos
    
	add_radiation_quantities!(eos, opa, compute_ross=compute_ross)
	fill_nan!(aos, opa)
	set_limits!(aos, opa, small=1e-30, large=1e30)

	eos, opa
end


function add_radiation_quantities!(eos, opa; compute_ross=false)
	S = similar(opa.κ)
	T = exp.(eos.lnT)
	for j in eachindex(eos.lnRho)
		S[:, j, :] .= TSO.Bλ(opa.λ, T)
	end

	if compute_ross
		rosseland_opacity!(eos.lnRoss, @axed(eos), opa)
		transfer_rosseland!(@axed(eos), opa)
	end

	eos, opa
end