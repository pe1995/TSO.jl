### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 89cb2c56-14bd-11ee-07ac-714f17e31c1c
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); Pkg.add("TimerOutputs")
	using TSO
	using ProfileVega
	using Profile
	using Plots
	using Optim
	using BenchmarkTools
	using TimerOutputs
end

# ╔═╡ d1741019-41c9-4701-85cb-96af9db2a4c0
md"# Optimizing Opacity Binning"

# ╔═╡ 3366cfb9-aa10-42e9-96a9-d1873038ac41
md"## Downsample opacity table
One major bottleneck is the size of the unbinned opacity table (104x104x100000). It should be possible however to downsample it by a factor of ~10 during the iterations of the optimization to find the ideal binning, and then use this binning scheme for the entire table. A first attempt would be to simply pick only every 10th point in the table."

# ╔═╡ 7f9748ef-c9dc-4063-b8b6-4f8adabb88aa
big_table = "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6"

# ╔═╡ 0fbc59d1-26e0-44d3-8a40-06c3efa536a7
eos = reload(SqEoS, joinpath(big_table, "combined_ross_eos.hdf5"))

# ╔═╡ 12a942b4-cf2a-40fb-bab4-efaad4062949
opacity = reload(
	SqOpacity, 
	joinpath(big_table, "combined_opacities.hdf5"), 
	mmap=true
)

# ╔═╡ 11ce2612-0f83-46b9-96f7-03cebff165d6
lopacity = reload(
	SqOpacity, 
	joinpath(big_table, "combined_Lopacities.hdf5"), 
	mmap=true
)

# ╔═╡ 3f0d7546-6c01-4b5a-a576-67de888e0a97
copacity = reload(
	SqOpacity, 
	joinpath(big_table, "combined_Copacities.hdf5"), 
	mmap=true
)

# ╔═╡ a128aace-a5e9-418d-ab08-e958368d5068
fopacity = reload(
	SqOpacity, 
	joinpath(big_table, "combined_formation_opacities.hdf5"), 
	mmap=true
)

# ╔═╡ ff134b65-7801-41a7-afcf-36b66bffd8db
# ╠═╡ show_logs = false
model = @optical Average3D(eos, "sun_stagger.dat") eos opacity

# ╔═╡ e2333b49-2071-43fe-aa27-6646e8403781
"""
	downsample_opacity(table; pick_every=10)

Downsample the opacity `table` by picking every `pick_every` point in the opacity space.
"""
function downsample_mask(line, cont, eos, model; pick_cont=20, pick_line=2)	
	kline = lookup(eos, line, :κ, model.lnρ, model.lnEi)
	kcont = lookup(eos, cont, :κ, model.lnρ, model.lnEi)

	i0 = log.(model.τ) .< 1.0
	mask = falses(length(line.λ))	

	j = 1
	k = 1
	for i in eachindex(mask)
		mask[i] = if TSO.median(kline[i0, i]) .> TSO.median(kcont[i0, i])
			true
		else
			false
		end
		
		mask[i] = if mask[i]
			if k == pick_line
				k=1
				true
			else
				k+=1
				false
			end
		else
			if j == pick_cont
				j=1
				true
			else
				j+=1
				false
			end
		end
	end

	
	mask
end

# ╔═╡ 7114478a-5a22-4bbb-9c77-e3b4c8de8de4
m = downsample_mask(lopacity, copacity, eos, model, pick_cont=100, pick_line=30)

# ╔═╡ c904317e-580a-4512-ad66-ff54ffd0c6b6
count(m)

# ╔═╡ 28aa0a5f-2f81-4a70-abcf-e823e0c5a176
"""
	downsample_opacity(table; pick_every=10)

Downsample the opacity `table` by picking every `pick_every` point in the opacity space.
"""
function downsample_opacity(table; pick_every=10, m=nothing)	
	m = if isnothing(m) 
		1:pick_every:length(table.λ)
	else
		m
	end
	λ   = table.λ[m]
	κ   = table.κ[:, :, m]
	src = table.src[:, :, m]
	κ_ross = deepcopy(table.κ_ross)
	
	SqOpacity(κ, κ_ross, src, λ, table.optical_depth)
end

# ╔═╡ e7635ae8-015b-4dc3-a2bf-a5b5f2479b0d
"""
	downsample_opacity(table; pick_every=10)

Downsample the opacity `table` by picking every `pick_every` point in the opacity space.
"""
function downsample_formation_opacity(table; pick_every=10, m=nothing)	
	m = if isnothing(m) 
		1:pick_every:length(table.λ)
	else
		m
	end
	λ   = table.λ[m]
	κ   = table.κ[m]
	src = table.src[:, :, m]
	κ_ross = table.κ_ross[m]
	
	SqOpacity(κ, κ_ross, src, λ, table.optical_depth)
end

# ╔═╡ 99126871-98ab-4ce6-bc62-4dc9402c40d7
opacity_small = downsample_opacity(opacity, m=m)

# ╔═╡ e234ad2c-e223-47d7-af80-f37b3471480c
fopacity_small = downsample_formation_opacity(fopacity, m=m)

# ╔═╡ d59cc9a9-511e-45de-a877-a543e7ae0ec1
md"## Binning speed
One way to optimize the binning is to optimize the time it takes to bin opacities. For this we use one of the default bins"

# ╔═╡ b98716f6-6bbc-4942-9ad0-43873fd275f5
quadrants = [
	TSO.Quadrant((0.0, 4.0), (-100, 0.5), 1),
	TSO.Quadrant((0.0, 4.0), (0.5, 4.0), 1),
	TSO.Quadrant((0.0, 4.0), (4.0, 100.), 1),
	TSO.Quadrant((4.0, 100.0), (-100, 1.0), 1),
	TSO.Quadrant((4.0, 100.0), (1.0, 100), 1)
]

# ╔═╡ 793a4ccc-9cc9-4c29-8b10-278f18d0f275
bins = binning(
	ClusterBinning(
		TSO.KmeansBins;
	    opacities=opacity_small, 
	    formation_opacity=-log10.(fopacity_small.κ_ross), 
	    Nbins=5, 
		quadrants=quadrants
	),
	opacity_small,
	-log10.(fopacity_small.κ_ross)
)

# ╔═╡ 1bc86853-ccc5-4b12-9340-592b3e3163fc
begin
	plot(framestyle=:box, grid=false)
	
	scatter!(
		log10.(opacity_small.λ), 
		-log10.(fopacity_small.κ_ross), 
		marker_z=TSO.assignments(bins),
		ms=3,
		color=:rainbow, markeralpha=0.4, label=nothing
	)
end

# ╔═╡ 134556db-0b34-4164-a586-d0b68bde3713
md"The actual binning requires some time..."

# ╔═╡ 9ac68b7d-d4ec-4b2a-88ac-1f16abb68123
weights = ω_midpoint(opacity_small)

# ╔═╡ b999c153-1bc1-4d88-a6a9-a96924bf3271
# ╠═╡ show_logs = false
@benchmark binned_small = tabulate(
	bins, weights, eos, opacity_small, transition_model=model
)

# ╔═╡ a57402ef-9feb-400f-acf4-2e8dd3f7cd66
md"Lets see if we can make it faster!"

# ╔═╡ ac43d5dd-968c-4c5e-8403-e833bd391593
begin
	function do_binning!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
						 ρ, λ, binning, Temp, weights, κ, κ_scat)
		rhoBins = size(χBox, 2)
		AxBins = size(χBox, 2)
	
		bnu  = Ref(Float32(0.0))
	    dbnu = Ref(Float32(0.0))
	    b ::Int = 0

		B  .= 0.0
	    δB .= 0.0
		
		@inbounds for k in eachindex(λ)
			b = binning[k]
			@inbounds for j in 1:rhoBins
		        @fastmath @inbounds for i in 1:AxBins
	                TSO.Bν!(bnu,   λ[k], Temp[i, j]) 
	                TSO.δBν!(dbnu, λ[k], Temp[i, j]) 
	
	                χBox[ i, j, b] += weights[k] * κ[i, j, k] * bnu[]  
	            	χRBox[i, j, b] += weights[k] * 1.0 / κ[i, j, k] * dbnu[]
	                
					χ_thin[i,j, b] = χRBox[i, j, b]
					
	                SBox[ i, j, b]  += weights[k] * bnu[] 
	                B[ i, j, b]     += weights[k] * bnu[]  
	                δB[i, j, b]     += weights[k] * dbnu[]
	                κBox[ i, j, b]  += weights[k] * weights[k] * 
												(κ[i, j, k] - κ_scat[i, j, k]) * bnu[]
	            end
	        end
	    end

		χBox   ./= B
	    κBox   ./= B
	    χRBox  ./= δB
	    χ_thin ./= δB

		@inbounds for j in 1:rhoBins
			χBox[:,  j, :] .*= ρ[j]
			κBox[:,  j, :] .*= ρ[j]
			χRBox[:, j, :] ./= ρ[j]
		end
	
	    χRBox[χRBox .<= 1e-30] .= 1e-30
	    χRBox[χRBox .>= 1e30]  .= 1e30
	    χBox[ χBox  .<= 1e-30] .= 1e-30
	    χBox[ χBox  .>= 1e30]  .= 1e30
	    SBox[SBox   .<= 1e-30] .= 1e-30
	    SBox[SBox   .>= 1e30]  .= 1e30
	end
	
	function do_binning!(B, δB, SBox, κBox, χBox, χRBox, χ_thin,
						 ρ, λ, binning, Temp, weights, κ)

		rhoBins = size(χBox, 2)
		AxBins = size(χBox, 2)
	
		bnu  = Ref(Float32(0.0))
	    dbnu = Ref(Float32(0.0))
	    b ::Int = 0

		B  .= 0.0
	    δB .= 0.0
		
		@inbounds for k in eachindex(λ)
			b = binning[k]
			@inbounds for j in 1:rhoBins
		        @fastmath @inbounds for i in 1:AxBins
	                TSO.Bν!(bnu,   λ[k], Temp[i, j]) 
	                TSO.δBν!(dbnu, λ[k], Temp[i, j]) 
	
	                χBox[ i, j, b] += weights[k] * κ[i, j, k] * bnu[]  
	            	χRBox[i, j, b] += weights[k] * 1.0 / κ[i, j, k] * dbnu[]
	                
					χ_thin[i,j, b] = χRBox[i, j, b]
					
	                SBox[ i, j, b]  += weights[k] * bnu[] 
	                B[ i, j, b]     += weights[k] * bnu[]  
	                δB[i, j, b]     += weights[k] * dbnu[]
	                κBox[ i, j, b]  += weights[k] * weights[k] * κ[i, j, k] * bnu[]
	            end
	        end
	    end

		χBox   ./= B
	    κBox   ./= B
	    χRBox  ./= δB
	    χ_thin ./= δB

		@inbounds for j in 1:rhoBins
			χBox[:,  j, :] .*= ρ[j]
			κBox[:,  j, :] .*= ρ[j]
			χRBox[:, j, :] ./= ρ[j]
		end
	
	    χRBox[χRBox .<= 1e-30] .= 1e-30
	    χRBox[χRBox .>= 1e30]  .= 1e30
	    χBox[ χBox  .<= 1e-30] .= 1e-30
	    χBox[ χBox  .>= 1e30]  .= 1e30
	    SBox[SBox   .<= 1e-30] .= 1e-30
	    SBox[SBox   .>= 1e30]  .= 1e30
	end
end

# ╔═╡ ee4a9907-62db-4797-9e66-caea923c5260
function bin_opacities(binning, weights, aos::E, opacities, scattering=nothing; 
						transition_model=nothing) where {E<:AxedEoS}
    eos = aos.eos
    eaxis = is_internal_energy(aos)

    iscat = isnothing(scattering)
    remove_from_thin = !iscat

    radBins = length(unique(binning))
    rhoBins = length(eos.lnRho)
    AxBins  = aos.energy_axes.length
    T       = eltype(eos.lnRho)

    χBox   = zeros(T, AxBins, rhoBins, radBins)    
    χRBox  = zeros(T, AxBins, rhoBins, radBins)
    χ_thin = similar(χRBox)
    κBox   = zeros(T, AxBins, rhoBins, radBins)
    SBox   = zeros(T, AxBins, rhoBins, radBins)
    Temp   = zeros(T, AxBins, rhoBins)
	
    if ndims(eos.lnT) == 2
        Temp .= exp.(eos.lnT)
    else
        TSO.puff_up!(Temp, exp.(eos.lnT))
    end

    B  = zeros(T, AxBins, rhoBins, radBins)
    δB = zeros(T, AxBins, rhoBins, radBins)

    bnu  ::Float32  = 0.0
    dbnu ::Float32  = 0.0
    b    ::Int      = 0

    ρ = exp.(eos.lnRho)

    # None of the bins are allowed to be empty as this stage
    @assert all(binning .!= 0)

	# Do stuff
	if (!iscat) 
		do_binning!(
			B, δB, SBox, κBox, χBox, χRBox, χ_thin, ρ,
			opacities.λ, binning, Temp, weights, 
			opacities.κ, scattering.κ_scat
		)
	else
		do_binning!(
			B, δB, SBox, κBox, χBox, χRBox, χ_thin, ρ,
			opacities.λ, binning, Temp, weights, 
			opacities.κ
		)
	end
	
	κ_ross = T(1.0) ./χRBox

    wthin, wthick = if isnothing(transition_model)
        wthin_l = exp.(T(-1.5e7) .* κ_ross)
        wthick_l = T(1.0) .- wthin_l
        wthin_l, wthick_l
    else
        pure_rosseland = TSO.RegularOpacityTable(
			κ_ross, opacities.κ_ross, SBox, collect(T, 1:radBins), false
		)
        TSO.opacity_transition(aos, pure_rosseland, transition_model)
    end

    ## Opacity
    opacity_table = remove_from_thin ? 
                    wthin .* κBox .+ wthick .* κ_ross : 
                    wthin .* χBox .+ wthick .* κ_ross

    ## ϵ table
    ϵ_table = κBox ./ χBox

    ## Source function table --> S                = x/x+s * B + s/x+s * J 
    ##                       --> thermal emission = x/x+s * B
    S_table = SBox

    TSO.BinnedOpacities(
		SqOpacity(
			opacity_table, opacities.κ_ross, S_table,
			collect(T, 1:radBins), false
		), 
		ϵ_table, 
		wthin
	)
end

# ╔═╡ 1a734578-7586-4340-92f2-879d1b285d88
begin
	#reset_timer!(to)
	
	binned_fast = bin_opacities(
		bins, weights, @axed(eos), opacity_small, transition_model=model
	)

	#show(to)
end

# ╔═╡ 0bc57b5d-79d1-456d-8e11-2380d9c82582
md"## Quantify goodness of binning"

# ╔═╡ 720ec7e7-0ce9-4f78-9308-47fe5bc5cea2
md"# Radiative transfer
Compute binned and unbinned radiative transfer and compare the results."

# ╔═╡ c9faaf88-b794-49e3-b4fe-9afc8c28282e
md"Applying that the opa opacities are in fact binned"

# ╔═╡ 81634146-c8a4-4d2b-b6ae-a4426ef6e2c3
opacities = binned_fast

# ╔═╡ 54b9f7be-f8be-4095-88ab-ccd3a0f6001d
md"Applying that the opa_raw opacities are in fact unbinned"

# ╔═╡ ed5b10cd-3b9f-4b30-b6ce-cc840adbb21e
opacities_raw = @binned opacity_small eos

# ╔═╡ af810105-2916-4054-9b8b-f138cbf621cc
weights_full = ω_midpoint(opacity_small)

# ╔═╡ 503eb45c-9ce4-421c-b703-e7b9766ce25d
solver_raw = Solver(model, @axed(eos), opacities=opacities_raw)

# ╔═╡ 55561895-edf9-43b0-852d-8a427aa6030c
solver = Solver(model, @axed(eos), opacities=opacities)

# ╔═╡ aea89619-4293-4430-b399-45eb2649ef64
md"Solve the radiative transfer using those solvers"

# ╔═╡ 762e5c41-0963-48b3-9a91-f910353b9f10
q = Qr(solver)

# ╔═╡ d98c2e3f-c92f-4d02-95ae-c12d2f19ee6a
md"Additionally pass the weigts to the raw solver, because the integration over λ is not perfomed yet in the unbinned case."

# ╔═╡ 933ff2ad-5959-4318-a030-5c25ce34efa0
q_raw = Qr(solver_raw, weights_full)

# ╔═╡ ba53555b-d1c1-474d-91ad-c2d6d2541bae
z, lnT, τ = solver_raw.model[:, 1], solver_raw.model[:, 2], reverse(model.τ)

# ╔═╡ a8718882-93bd-43e1-9191-793becb0e462
begin
	mask = log10.(τ) .< 5

	plot(framestyle=:box, grid=false)
	
	scatter!(
		log10.(τ[mask]), q_raw[mask], 
		label="unbinned", color=:cyan
	)

	plot!(
		log10.(τ[mask]), q[mask], 
		label="binned", color=:black, lw=3
	)
end

# ╔═╡ 3ff70c92-e6cb-4d07-b09f-2379bff99368
begin
	plot(framestyle=:box, grid=false)
	
	plot!(
		log10.(τ[mask]), (q[mask] - q_raw[mask]) ./q_raw[mask], 
		color=:black, lw=3
	)

	plot!(ylim=(-0.5,0.5))
end

# ╔═╡ 02be6ad5-e7a4-408b-8d2d-36242a811eb3
"""
	goodness(binning, unbinned_heating)

Evaluate the goodness of the given heating based on the result for the unbinned table.
"""
function goodness(binned_opacity, unbinned_heating; model, eos, τ_min=-5, τ_max=0)
	solver = Solver(model, @axed(eos), opacities=binned_opacity)
	q = Qr(solver)

	τ = reverse(model.τ)
	mask = τ_min .< log10.(τ) .< τ_max

	log(sum(((q[mask] - unbinned_heating[mask])).^2))
end

# ╔═╡ c562ff9f-e73e-4849-a03b-53431986ca8c
goodness(binned_fast, q_raw, model=model, eos=eos)

# ╔═╡ 48656440-46f9-4846-a208-f4dbf70bd806
quadrants

# ╔═╡ 8bdb30bb-a6d1-405e-9329-cacd53c2c560
begin
	qi = deepcopy(quadrants)

	ledge1 = 4.0
	kedge1 = 0.5
	kedge2 = 4.0
	kedge3 = 1.0

	edges = [ledge1, kedge1, kedge2, kedge3]
	lower = [3.0, 0.1, 3.0, 0.5]
	upper = [4.5, 2.0, 5.0, 1.5]
	
	assignments(e) = begin
		[TSO.Quadrant((0.0, e[1]), (-100, e[2]), 1),
		TSO.Quadrant((0.0, e[1]), (e[2], e[3]), 1),
		TSO.Quadrant((0.0, e[1]), (e[3], 100.), 1),
		TSO.Quadrant((e[1], 100.0), (-100, e[4]), 1),
		TSO.Quadrant((e[1], 100.0), (e[4], 100), 1)]
	end
end

# ╔═╡ a2dadf50-1885-44ee-b219-a2e116aa04e9
begin
	"""
	Fitting the edges of given Quadrants.
	Detect based on the Quadrants which edges need to be fitted.
	"""
	mutable struct QuadrantEdges
		quadrants
		edges
		assignments
	end

	set_edges!(q::QuadrantEdges, edges) = begin
		q.edges .= edges
		q.quadrants = assignments(q.edges)
	end
end

# ╔═╡ 0d22b2f0-f04b-4e96-8ab8-9c75e1cec519
"""
	likelihood(bin_edges, heating, model, eos, opacities; kwargs...)

From the given `bin_edges` compute the binning, bin the `opacities` and compare the results to the given `heating`.
"""
function likelihood(bin_edges, 
						heating, 
						model, 
						eos, 
						opacities, 
						formation_opacities, 
						quads; 
						kwargs...)
	# set the edges and compute the quadrants
	e = [ledge1, kedge1, kedge2, kedge3]
	de = bin_edges .- e

	if any((e .+ 100.0 .* de) .> upper) | any((e .+ 100.0 .* de) .< lower)
		return Inf
	end

	quadrants = quads(e .+ 100.0 .* de)
	
	#@show e .+ 50.0 .* de
	# compute the binning
	bins = binning(
		ClusterBinning(
			TSO.KmeansBins;
		    opacities=opacities, 
		    formation_opacity=-log10.(formation_opacities.κ_ross), 
		    Nbins=5, 
			quadrants=quadrants
		),
		opacities,
		-log10.(formation_opacities.κ_ross)
	)

	mask = bins .== 0
	if any(mask) 
		@warn "some bins are 0 ($(edges)), N=$(count(mask))"
		
		for i in eachindex(mask)
			if mask[i]
				j=i
				while mask[j]
					j-=1
				end
				bins[i] = bins[j]
			end
		end
	end

	# integration weights
	weights = ω_midpoint(opacities)

	# bin the opacities
	binned_opacity = bin_opacities(
		bins, weights, @axed(eos), opacities, transition_model=model
	)

	# estimate the goodness of this binning (log)
	goodness(binned_opacity, heating, eos=eos, model=model, kwargs...)
end

# ╔═╡ 1d9c608a-73d5-4a31-85dd-659bf5690c07
llhood(edges) = likelihood(
	edges, q_raw, model, eos, opacity_small, fopacity_small, assignments
) 

# ╔═╡ e76b5a26-34cf-402c-a7bf-492fa704855c
llhood(edges)

# ╔═╡ 68a4e23b-770e-487a-84a2-d3e6259b7493
begin
	p = optimize(
		llhood, 
		lower, upper, deepcopy(edges), 
		Fminbox(GradientDescent()), 
		Optim.Options(allow_f_increases=true, iterations=10)
	)
	
	@show p
	best_fit = Optim.minimizer(p)
end

# ╔═╡ b5fa60c0-bef1-4c20-b292-6c2d41ab5657
#best_fit

# ╔═╡ Cell order:
# ╟─d1741019-41c9-4701-85cb-96af9db2a4c0
# ╠═89cb2c56-14bd-11ee-07ac-714f17e31c1c
# ╟─3366cfb9-aa10-42e9-96a9-d1873038ac41
# ╠═7f9748ef-c9dc-4063-b8b6-4f8adabb88aa
# ╠═0fbc59d1-26e0-44d3-8a40-06c3efa536a7
# ╠═12a942b4-cf2a-40fb-bab4-efaad4062949
# ╠═11ce2612-0f83-46b9-96f7-03cebff165d6
# ╠═3f0d7546-6c01-4b5a-a576-67de888e0a97
# ╠═a128aace-a5e9-418d-ab08-e958368d5068
# ╠═ff134b65-7801-41a7-afcf-36b66bffd8db
# ╟─e2333b49-2071-43fe-aa27-6646e8403781
# ╠═7114478a-5a22-4bbb-9c77-e3b4c8de8de4
# ╠═c904317e-580a-4512-ad66-ff54ffd0c6b6
# ╠═28aa0a5f-2f81-4a70-abcf-e823e0c5a176
# ╠═e7635ae8-015b-4dc3-a2bf-a5b5f2479b0d
# ╠═99126871-98ab-4ce6-bc62-4dc9402c40d7
# ╠═e234ad2c-e223-47d7-af80-f37b3471480c
# ╟─d59cc9a9-511e-45de-a877-a543e7ae0ec1
# ╠═b98716f6-6bbc-4942-9ad0-43873fd275f5
# ╠═793a4ccc-9cc9-4c29-8b10-278f18d0f275
# ╟─1bc86853-ccc5-4b12-9340-592b3e3163fc
# ╟─134556db-0b34-4164-a586-d0b68bde3713
# ╠═9ac68b7d-d4ec-4b2a-88ac-1f16abb68123
# ╠═b999c153-1bc1-4d88-a6a9-a96924bf3271
# ╟─a57402ef-9feb-400f-acf4-2e8dd3f7cd66
# ╟─ac43d5dd-968c-4c5e-8403-e833bd391593
# ╟─ee4a9907-62db-4797-9e66-caea923c5260
# ╠═1a734578-7586-4340-92f2-879d1b285d88
# ╟─0bc57b5d-79d1-456d-8e11-2380d9c82582
# ╟─720ec7e7-0ce9-4f78-9308-47fe5bc5cea2
# ╟─c9faaf88-b794-49e3-b4fe-9afc8c28282e
# ╠═81634146-c8a4-4d2b-b6ae-a4426ef6e2c3
# ╟─54b9f7be-f8be-4095-88ab-ccd3a0f6001d
# ╠═ed5b10cd-3b9f-4b30-b6ce-cc840adbb21e
# ╠═af810105-2916-4054-9b8b-f138cbf621cc
# ╠═503eb45c-9ce4-421c-b703-e7b9766ce25d
# ╠═55561895-edf9-43b0-852d-8a427aa6030c
# ╟─aea89619-4293-4430-b399-45eb2649ef64
# ╠═762e5c41-0963-48b3-9a91-f910353b9f10
# ╟─d98c2e3f-c92f-4d02-95ae-c12d2f19ee6a
# ╠═933ff2ad-5959-4318-a030-5c25ce34efa0
# ╠═ba53555b-d1c1-474d-91ad-c2d6d2541bae
# ╟─a8718882-93bd-43e1-9191-793becb0e462
# ╟─3ff70c92-e6cb-4d07-b09f-2379bff99368
# ╟─02be6ad5-e7a4-408b-8d2d-36242a811eb3
# ╠═c562ff9f-e73e-4849-a03b-53431986ca8c
# ╠═a2dadf50-1885-44ee-b219-a2e116aa04e9
# ╠═48656440-46f9-4846-a208-f4dbf70bd806
# ╠═8bdb30bb-a6d1-405e-9329-cacd53c2c560
# ╠═0d22b2f0-f04b-4e96-8ab8-9c75e1cec519
# ╠═1d9c608a-73d5-4a31-85dd-659bf5690c07
# ╠═e76b5a26-34cf-402c-a7bf-492fa704855c
# ╠═68a4e23b-770e-487a-84a2-d3e6259b7493
# ╠═b5fa60c0-bef1-4c20-b292-6c2d41ab5657
