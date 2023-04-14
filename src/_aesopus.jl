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

    lnRho = Base.convert.(type, first(content) .* log(10))
    for it in 1:nT
        lnT[it] = Base.convert.(type, first(content[it+1]) .* log(10))

        for jr in 1:nR
            lnEi[it,   jr] = Base.convert(type, content[ it+1][jr+1] * log(10))
            lnRoss[it, jr] = Base.convert(type, contentK[it+1][jr+1] * log(10))
            lnNe[it,   jr] = Base.convert(type, contentN[it+1][jr+1] * log(10))
            lnPg[it,   jr] = Base.convert(type, contentP[it+1][jr+1] * log(10))
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

    ## Sort lambda from low to high
    reverse!(l)
    reverse!(κ, dims=3)

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




#= Utility functions =#

"""
Read and parse content of ascii files
"""
get_content(path, type) = open(path, "r") do f
    content = readlines(f)
    [parse.(type, strip.(split(strip(c)))) for c in content if length(split(strip(c))) > 1]
end
