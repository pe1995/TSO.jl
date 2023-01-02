#=================================================#
#======= AESOPUS 2.0 Opacity + EoS ===============#
#=================================================#

#= Constructors for the EoS =#

"""
    AesopusEoS(path)

Read an Aesopus EoS from a HDF5 file. 
Reads it directly into the squaregas format.

# Examples
```jldoctest
julia> AesopusEoS(path)
 RegularEoSTable(...)
```
"""
AesopusEoS(path) = AesopusEoS_hdf5(path)

"""
    AesopusEoS(; [type=Float32], paths...)

Read an Aesopus EoS from ascii files. 
Reads it directly into the squaregas format.

# Examples
```jldoctest
julia> AesopusEoS(energy=path1, kross=path2, ne=path3, pg=path4)
 RegularEoSTable(...)
```
"""
AesopusEoS(; type=Float32, paths...) = AesopusEoS_ascii(; type, paths...)

"""
Read table in HDF5 format.
"""
function AesopusEoS_hdf5(path)
    f = HDF5.h5open(path, "r")

    lnT    = log.(HDF5.read(f["T"]))
    lnEi   = HDF5.read(f["lnEi"])
    lnRoss = HDF5.read(f["lnKross"])
    lnNe   = HDF5.read(f["lnNe"])
    lnPg   = HDF5.read(f["lnPg"])
    lnRho  = HDF5.read(f["lnRho"])

    close(f)

    @assert all(size(lnPg) .== (length(lnT), length(lnRho)))

    RegularEoSTable(lnRho, lnT, lnEi, lnPg, lnRoss, lnNe)
end

"""
Read tables in ascii format.
"""
function AesopusEoS_ascii(; energy, kross, ne, pg, type=Float32)
    content = get_content(energy, type)

    ## The number of lines is equal to the number of T points +1 
    nT = length(content) - 1 
    nR = length(first(content))

    lnT    = zeros(type, nT)
    lnRho  = zeros(type,     nR)
    lnEi   = zeros(type, nT, nR)
    lnRoss = zeros(type, nT, nR)
    lnNe   = zeros(type, nT, nR)
    lnPg   = zeros(type, nT, nR)

    contentK = get_content(kross, type)
    contentN = get_content(ne, type)
    contentP = get_content(pg, type)

    lnRho = log.(first(content))
    for it in 1:nT
        lnT[it] = log(first(content[it+1]))

        for jr in 1:nR
            lnEi[it,   jr] = content[ it+1][jr+1]
            lnRoss[it, jr] = contentK[it+1][jr+1]
            lnNe[it,   jr] = contentN[it+1][jr+1]
            lnPg[it,   jr] = contentP[it+1][jr+1]
        end
    end

    @assert all(size(lnPg) .== (length(lnT), length(lnRho)))

    RegularEoSTable(lnRho, reverse(lnT), reverse(lnEi, dims=1), reverse(lnPg, dims=1), reverse(lnRoss, dims=1), reverse(lnNe, dims=1))
end




#= Constructors for the Opacity table =#

"""
    AesopusOpacity(eos; sigma=path_to_file, λ_nodes=5776)

Read the Aesopus opacity file in binary format.
Contains the absolute opacity (line + cont + scattering) in cm2/g
Format: log10(T) log10(rho) log10(lambda/Angstrom) log10(sigma)
"""
AesopusOpacity(eos; paths...) = AesopusOpacity_binary(eos; paths...)

"""
Read table in binary format. Ensure consistency with the EoS.
"""
function AesopusOpacity_binary(eos::RegularEoSTable; sigma, λ_nodes=5776)
    lt = length(eos.lnT)
    lr = length(eos.lnRho)
    ll = λ_nodes
    ltot = lt * lr * ll

    T = zeros(Float64, lt)
    r = zeros(Float64, lr)
    l = zeros(Float64, ll)
    κ = zeros(Float64, lt, lr, ll)
    line = zeros(Float64, 4)

    f = open(sigma, "r")

    for k in 1:lt
        for j in 1:lr
            for i in 1:ll
                T[k], r[j], l[i], κ[k, j, i] = read!(f, line)
            end
        end
    end

    close(f)

    ## log10 -> ln
    T .*= log(10)
    r .*= log(10)
    l .*= log(10)
    κ .*= log(10)

    ## Sort T from low to high
    reverse!(T)
    reverse!(κ, dims=1)

    ## Convert to the requested type
    if eltype(T) != eltype(eos.lnT)
        tp   = eltype(eos.lnT)
        Tnew = zeros(tp, size(T)...)
        rnew = zeros(tp, size(r)...)
        lnew = zeros(tp, size(l)...)
        κnew = zeros(tp, size(κ)...)
        
        @inbounds for k in eachindex(l)
            @inbounds for j in 1:lr
                @inbounds for i in 1:lt
                    κnew[i, j, k] = Base.convert(tp, κ[i, j, k])
                    lnew[k]       = Base.convert(tp, l[k])
                    Tnew[i]       = Base.convert(tp, T[i])
                    rnew[j]       = Base.convert(tp, r[j])
                end
            end
        end
    else
        Tnew = T
        rnew = r
        lnew = l
        κnew = κ
    end
        

    ## The grids have to match (for now, could be interpolated if needed)
    @assert all(Tnew .≈ eos.lnT)
    @assert all(rnew .≈ eos.lnRho)

    ## Compute the source function (not integrated)
    src = similar(κnew)
    @inbounds for k in eachindex(lnew)
        @inbounds for j in 1:lr
            @inbounds for i in 1:lt
                src[i, j, k] = Bν(exp(lnew[k]), exp(Tnew[i]))
            end
        end
    end

    ## Create the opacity table
    RegularOpacityTable(exp.(κnew), exp.(eos.lnRoss), src, exp.(lnew), false)
end




#= EoS options =#

"""
    complement(opacities, eos_old, eos_new; lnRoss=:opacity)

Complement a RegularOpacityTable with an Regular EoS Table. The Opacities will be interpolated
to the EoS grid. Rosseland opacities can either be interpolated or recomputed (if lnRoss == :recompute).
Idea:
    New eos has rho_new + E_new axis (where E_new can be either T or Ei)
    old eso has rho + E axis (see above)

    1.) Old: T grid, New: E grid

        We want to know the opacity, which is tablulated on rho,T
        for a given rho_new,Ei_new. Since the old table is already 
        on absolute physical quantities, we can look up directly using 
        rho_new[j], lnT_new[i, j], lnT_new is always 2D.

    2.) Old: E grid, New: E grid

        Both are tabulated on E grid. Since E scales might be different. To use the 
        lookup interface, we have to convert the new E axis to the old one, then lookup.
        Because the old one is tablulated on E, this requires a bisect search using the 
        rho_new and lnT_new, lnT_new is always 2D.

    3.) Old: E grid, New: T grid
        The same as 2.) needs to be done, however, the lnT_new will be 1D.

    4.) Old: T grid, New: T grid
        The same as 1.) needs to be done, however, the lnT_new will be 1D.

"""
function complement(opa::O, eos_old::E1, eos_new::E2; lnRoss=:opacity) where {O<:RegularOpacityTable, E1<:RegularEoSTable, E2<:RegularEoSTable}
    eaxis_old = EnergyAxis(eos_old)
    eaxis_new = EnergyAxis(eos_new)

    κ  = zeros(eltype(opa.κ), size(eos_new.lnPg)..., length(opa.λ))
    S  = zeros(eltype(opa.κ), size(eos_new.lnPg)..., length(opa.λ))
    κR = zeros(eltype(opa.κ), size(eos_new.lnPg)...)

    old_E    = eaxis_old.name == :lnEi   ## Is the axis of the old table Ei?
    lnTNew2D = ndims(eos_new.lnT) == 2   ## Is the temperature_new 2D (e.g. Energy grid)

    EOldForNewEOS  = zeros(eltype(eos_new.lnEi), size(eos_new.lnPg)...)
    emin,emax,_,_  = limits(eos_old)

    @info "Old Grid:", eaxis_old.name
    @info "New Grid:", eaxis_new.name

    # We get the lookup functions for the old table (ross)
    κR_old = if lnRoss == :opacity
        lookup_function(eos_old, opa, :κ_ross)
    elseif lnRoss == :eos_old
        lookup_function(eos_old, :lnRoss)
    else
        (args...) -> 0.0
    end


    for j in eachindex(eos_new.lnRho)
        for i in eachindex(eaxis_new.values)
            ## We need to know the energy or the temperature 
            ## this point corresponds to on the old table
            lnT = lnTNew2D ? eos_new.lnT[i, j] : eos_new.lnT[i]

            EOldForNewEOS[i, j] = if old_E
                bisect(eos_old, lnRho=eos_new.lnRho[j], lnT=lnT, lnEi=[emin, emax])
            else
                lnT
            end

            ## Next we need to lookup the stuff we want on this new grid
            κR[i, j] = κR_old(EOldForNewEOS[i, j], eos_new.lnRho[j])
        end
    end

    ## Now the wavelength dependent stuff
    @inbounds for k in eachindex(opa.λ)
        # We get the lookup functions for the old table
        κ_old = lookup_function(eos_old, opa, :κ, k)
        S_old = lookup_function(eos_old, opa, :src, k)

        @inbounds for j in eachindex(eos_new.lnRho)
            @inbounds for i in eachindex(eaxis_new.values)
                κ[i, j, k] = κ_old(EOldForNewEOS[i, j], eos_new.lnRho[j])
                S[i, j, k] = S_old(EOldForNewEOS[i, j], eos_new.lnRho[j])
            end
        end
    end


    if lnRoss == :eos_new
        κR .= eos_new.lnRoss
    elseif lnRoss == :eos_old
        κR .= exp.(κR)
    end

    RegularOpacityTable(κ, κR, S, opa.λ, opa.optical_depth)
end



#= Utility functions =#

"""
Read and parse content of ascii files
"""
get_content(path, type) = open(path, "r") do f
    content = readlines(f)
    [parse.(type, strip.(split(strip(c)))) for c in content if length(split(strip(c))) > 1]
end
