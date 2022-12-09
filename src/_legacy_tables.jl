#= 
Reading legacy Stagger EoS tables and convert them to the SquareGas format required by DISPATCH 
=#

mutable struct LegacyTable{I<:Integer, F<:AbstractFloat}
    Nvtab   ::I
    Nrtab   ::I
    Ma      ::I
    Na      ::I
    Nvar    ::I
    Nelem   ::I
    Nxopt   ::I
    Nxre    ::I
    date    ::String
    time    ::String
    lnrmin  ::F
    lnrmax  ::F
    lnvmin  ::F
    lnvmax  ::F
    ul     ::F
    ur     ::F
    ut     ::F
    iel     ::Vector{String}
    abund   ::Vector{F}
    ar      ::Matrix{F}
    tab     ::Array{F, 4}
    parameters::Vector{String}
end



#= Constructors =#
"""
    LegacyOpacities(path)

Read a legacy Stagger Opacity table and convert it to the TSO format. 

# Examples
```jldoctest
julia> opacites = LegacyOpacities("kaprhoT.tab")
 TSO.LegacyTable{Int32, Float32}(...)
````
"""
function LegacyOpacities(path)
    @assert isfile(path)
    
    ## Read the file step by step and fill in the arrays
    f = FortranFiles.FortranFile(path, convert="big-endian")

    ## dimensions of the arrays
    NTtab, Nrtab, Ma, Na, Nlam, Nelem, Nxopt, Nxre = read(f, (Int32, 8))
    
    ## Read the header
    date, time, lnrmin, lnrmax, lntmin, lntmax, rul, rur, rut, iel, abund, ar = read(f, 
                                                                                    FString{8}, 
                                                                                    FString{8}, 
                                                                                    Float32,
                                                                                    Float32,
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    (FString{4}, Nelem), 
                                                                                    (Float32, Nelem), 
                                                                                    (Float32, Ma, Na))
    ## Read the rest of the table
    tab = zeros(Float32, NTtab, Nrtab, 3, Nlam+1)
    for j=1:Nlam+1
        tab[:, :, :, j] .= read(f, (Float32, NTtab, Nrtab, 3))
    end

    close(f)

    l = LegacyTable(NTtab ,
                Nrtab ,
                Ma    ,
                Na    ,
                Nlam  ,
                Nelem ,
                Nxopt ,
                Nxre  ,
                trimstring(date)  ,
                trimstring(time)  ,
                lnrmin,
                lnrmax,
                lntmin,
                lntmax,
                rul   ,
                rur   ,
                rut   ,
                trimstring.(iel)   ,
                abund ,
                ar    ,
                tab,
                ["opacity [cm-1]", "thermal emission", "destruction probability"])

    t, r = grid(l)

    for i in eachindex(r)
        l.tab[:, i, 1, :] .+= r[i]
    end

    l
end


"""
    LegacyEoS(path)

Read a legacy Stagger Opacity table and convert it to the TSO format. 

# Examples
```jldoctest
julia> opacites = LegacyEoS("kaprhoT.tab")
 TSO.LegacyTable{Int32, Float32}(...)
````
"""
function LegacyEoS(path)
    @assert isfile(path)
    
    ## Read the file step by step and fill in the arrays
    f = FortranFiles.FortranFile(path, convert="big-endian")

    Netab, Nrtab, Ma, Na, Ntbvar, Nelem, Nxopt, Nxrev = read(f, (Int32, 8))

    ## Read the header
    date, time, lnrmin, lnrmax, lnemin, lnemax, ul, ur, ut, iel, abund, ar = read(f, 
                                                                                    FString{8}, 
                                                                                    FString{8}, 
                                                                                    Float32,
                                                                                    Float32,
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    Float32, 
                                                                                    (FString{4}, Nelem), 
                                                                                    (Float32, Nelem), 
                                                                                    (Float32, Ma, Na))
    ## Read the rest of the table
    tab = read(f, (Float32, Netab, Nrtab, 3, Ntbvar))
    
    close(f)

    LegacyTable(Netab ,
                Nrtab ,
                Ma    ,
                Na    ,
                Ntbvar  ,
                Nelem ,
                Nxopt ,
                Nxrev  ,
                trimstring(date)  ,
                trimstring(time)  ,
                lnrmin,
                lnrmax,
                lnemin,
                lnemax,
                ul   ,
                ur   ,
                ut   ,
                trimstring.(iel)   ,
                abund ,
                ar    ,
                tab,
                ["lnPtot", "lnkap_Ross", "lnT", "lnkap_Planck","lnkap_5000"])    
end



## Make grids from min/max of axes
grid(t::LegacyTable) = begin
    r = range(t.lnrmin, t.lnrmax, length=t.Nrtab) |> collect
    v = range(t.lnvmin, t.lnvmax, length=t.Nvtab) |> collect

    v,r
end

grid(f::Function, t::LegacyTable) = begin
    a, b = grid(t)
    f.(a), f.(b)
end



## Conversion to the TSO API
##  - Note that this includes interpolation of -
##  - the binned opacity table to the E grid!  -
"""
    toTSO(eos, opacites)

Convert a legacy EoS + opacity table to the TSO API, so that it can 
be ported to DISPATCH easily. This involves interpolation of the 
opacities to the EoS grid, since in legacy tables opacities and source
function are listed on a temperature grid. Howerver in SquareGas table
a regular grid in energy is required! This may include extrapolation.

Please also note that this function assumes that: \n
    EoS:
        1.) 4th entry is electron density
        2.) 1st entry is gas pressure. Note that the doc says this 
            is total pressure, including radiation! For rho=1e-7 
            the differences compare to TSO where small below logT=10.5.
            However, this might require further checking!
        3.) The Rosseland optical depth from the EoS is used (supposed to be cm-1)
    Opacity:
        1.) First entry is opacity
        2.) Second entry is source function
        3.) It is valid to perform linear extrapolation in temperature 
            to cover the EoS.

# Examples
```jldoctest
julia> toTSO(eos, opacites)
 TSO.RegularEoSTable(...), TSO.RegularOpacityTable(...)
```
"""
function toTSO(eos::LegacyTable, opacities::LegacyTable)
    ## new axes
    newE, newR = grid(eos)
    oldT, oldR = grid(opacities)

    T = eltype(opacities.tab)

    nbins = size(opacities.tab, 4) - 1
    λ     = collect(T, 1:nbins)
    κ = zeros(T, length(newE), length(newR), nbins)
    S = similar(κ)
    ϵ = similar(S)
    ϵ .= 0.0
    RossCM2G = deepcopy(eos.tab[:, :, 1, 2])

    for i in eachindex(newR)
        T = view(eos.tab, :, i, 1, 3)
        for bin in 1:nbins
            κ[:, i, bin] .= legacyLookup(opacities, 1, newR[i], T, bin)
            S[:, i, bin] .= legacyLookup(opacities, 2, newR[i], T, bin)
        end
        RossCM2G[:, i] .-= newR[i]
    end

    RegularEoSTable(newR, eos.tab[:, :, 1, 3], newE, eos.tab[:, :, 1, 1], RossCM2G, eos.tab[:, :, 1, 4]),
    RegularOpacityTable(κ, ϵ, S, λ, false)
end

function legacyLookup(opacities, what, rho, T, args...)
    Tg, rg = grid(opacities)
    ip = extrapolate(interpolate((Tg, rg), view(opacities.tab, :, :, what, args...), Gridded(Linear())), Line())
    ip.(T, rho)
end
