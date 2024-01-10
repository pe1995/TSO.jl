#= Running M3D =#

function opacityTableInput(modelatmosfolder; 
                            lnT = range(log(1.1e3), log(5.5e5); length=59) |> collect, 
                            lnρ = range(log(1e-15), log(1e-3); length=59) |> collect)
	!isdir(modelatmosfolder) && mkdir(modelatmosfolder)

	z = range(-1, 1, length=length(lnT)) |> collect
	v = similar(lnT)

	paths = String[]
	for i in eachindex(lnρ)
		path_new = joinpath(modelatmosfolder, "TSO-M3D-$i")
		saveAsText(path_new, T=exp.(lnT), ρ=fill!(v, exp.(lnρ[i])), z=z)

		append!(paths, ["TSO-M3D-$i"])
	end

	paths
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
	c = run.run.read_patch_save("chi", concat=false)
	s = run.run.read_patch_save("src", concat=false)
	e = run.run.read_patch_save("E", concat=false)
	p = run.run.read_patch_save("Pg", concat=false)
	n = run.run.read_patch_save("Ne", concat=false)
	
	ck = pyconvert(Array, c.keys())
	sk = pyconvert(Array, s.keys())
	ek = pyconvert(Array, e.keys())
	pk = pyconvert(Array, p.keys())
	nk = pyconvert(Array, n.keys())
	
	chi = pyconvert(Array, c[ck[1]])[:, :, :, :, 1]
	src = pyconvert(Array, s[sk[1]])[:, :, :, :, 1]
	E = pyconvert(Array, e[ek[1]])[:, :, :, 1, 1]
	P = pyconvert(Array, p[pk[1]])[:, :, :, 1, 1]
	N = pyconvert(Array, n[nk[1]])[:, :, :, 1, 1]

	r = run.rho
	T = run.temp
	l = run.lam

	chi[chi .< 1e-30] .= NaN
	chi[chi .> 1e30] .=  NaN
	src[src .< 1e-30] .= NaN
	src[src .> 1e30] .= NaN
	
	E[E .< 1e-30] .= NaN
	E[E .> 1e30] .= NaN

	P[P .< 1e-30] .= NaN
	P[P .> 1e30] .= NaN

	N[N .< 1e-30] .= NaN
	N[N .> 1e30] .= NaN
	
	l, r, T, E, P, N, chi, src
end

"""
    collect_opacity(runs)

Extract EoS + Opacities from given runs. Collect them together in a table.
"""
function collect_opacity(runs)
	l, r, T, E, P, N, chi, src = [], [], [], [], [], [], [], []
	for i in eachindex(runs)
		li, ri, Ti, Ei, Pi, Ni, chii, srci = get_opacity(runs[i])

		append!(l, [li])
		append!(r, [ri])
		append!(T, [Ti])
		append!(E, [Ei])
		append!(P, [Pi])
		append!(N, [Ni])
		append!(chi, [chii])
		append!(src, [srci])
	end

	# We sort everything into arrays (TxRho)
	Rho = Iterators.flatten(r) |> collect |> unique |> sort
	Temp = Iterators.flatten(T) |> collect |> unique |> sort
	lam = first(l)

	@assert length(Rho)==length(runs)
	@assert length(Temp)==length(first(T))
	@assert Temp == first(T)

	te = eltype(first(E))
	to = eltype(first(chi))

	lam = to[lam...]
	
	lnEi = zeros(te, length(Temp), length(Rho))
	lnPg = zeros(te, length(Temp), length(Rho))
	lnNe = zeros(te, length(Temp), length(Rho))
	lnChi = zeros(te, length(Temp), length(Rho), length(lam))
	lnSrc = zeros(te, length(Temp), length(Rho), length(lam))

	for i in eachindex(r)
		ri = first(r[i])
		ir = findfirst(x->x≈ri, Rho)
		
		lnEi[:, ir] .= log.(E[i][1, 1, :])
		lnPg[:, ir] .= log.(P[i][1, 1, :])
		lnNe[:, ir] .= log.(N[i][1, 1, :])
		
		lnChi[:, ir, :] .= log.(chi[i][1, 1, :, :])
		#lnSrc[:, ir, :] .= log.(src[i][1, 1, :, :])
		lnSrc[:, ir, :] .= log.(Base.convert.(to, Bλ(lam, Temp)))
	end

	eos = TSO.SqEoS(
		log.(Rho), log.(Temp), lnEi, lnPg, zeros(to, length(Temp), length(Rho)), lnNe
	)

	opa = TSO.SqOpacity(
		exp.(lnChi), ones(to, length(Temp), length(Rho)), exp.(lnSrc), lam, false
	)

    aos = @axed eos

    fill_nan!(aos, opa)
    set_limits!(aos, opa)

	eos, opa
end