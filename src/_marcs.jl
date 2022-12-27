#================================================
========== MARCS Sampled Opacities ==============
================================================#

#= Intermediate opacity table type =#

"""
    MARCSOS(T, pe, pg, ρ, κ_cont, κ_scatter, κ_line_atom, κ_line_mol)

Content of a MARCS Opacity sampling file. Intermediate step before converting
to TS format.
"""
struct MARCSOS{F<:AbstractFloat,I<:Integer,
                NT,NPe,NPg,NR, NK}
    T         ::Array{F,NT}
    pe        ::Array{F,NPe}
    pg        ::Array{F,NPg}
    ρ         ::Array{F,NR}
    κ_c       ::Array{F,NK}
    κ_s       ::Array{F,NK}
    κ_la      ::Array{F,NK}
    κ_lm      ::Array{F,NK}
    λ         ::Vector{F}
    ω         ::Vector{F}
    abundance ::Dict{I, F}
end




#= Constructors for the Opacity tables =#

"""
    MARCSOpacity(paths...)

Read OS from MARCS (prepared by Bertrand Plez) in legacy format.
The relevant units are: \n
    integer nlambda, ntau, j, k, i, natoms
    real weight(150000),lambda(150000)
    real abscont(150),abslineat(150),scattcont(150)
    real abslinemol(150)
    real t(150),pe(150),pg(150),rho(150)
    integer species(100)
    real abundlog(100)
    character header_text*60,input_file*128
Read all the tables in an intermediate format.
"""
function MARCSOpacity(paths::String...)
    mos = []
    for path in paths
        append!(mos, [MARCSOS(path)])
    end
    
    mos
end

"""
    MARCSOS(path)

Read a single MARCS OS file and save it as an intermediate type.
"""
function MARCSOS(path::String)
    f = FortranFiles.FortranFile(path, "r")
        
    ## Read the header
    header = read(f, FortranFiles.FString{60}) |> trimstring
    while (header[1:1]=="*") & (header[1:2] != "**")
        header = read(f, FortranFiles.FString{60}) |> trimstring
        #println(header)
    end

    ## Read the data
    natoms  = read(f, Int32)
    species = read(f, (Int32, natoms))
    abunds  = read(f, (Float32, natoms))
    nlambda, ntau = read(f, Int32, Int32)

    T  = zeros(Float32, ntau)
    pe = zeros(Float32, ntau)
    pg = zeros(Float32, ntau)
    ρ  = zeros(Float32, ntau)
  
    eos = read(f, (Float32, 4, ntau))
    T, pe, pg, ρ = eos[1, :],  eos[2, :],  eos[3, :],  eos[4, :]

    κ_c  = zeros(Float32, ntau, nlambda)
    κ_s  = zeros(Float32, ntau, nlambda)
    κ_la = zeros(Float32, ntau, nlambda)
    κ_lm = zeros(Float32, ntau, nlambda)
    λ    = zeros(Float32, nlambda)
    ω    = zeros(Float32, nlambda)
    for j in 1:nlambda
        λ[j], ω[j], κ = read(f, Float32, Float32, (Float32, 4, ntau))
        κ_c[:,  j] = κ[1, :] 
        κ_s[:,  j] = κ[2, :]
        κ_la[:, j] = κ[3, :]
        κ_lm[:, j] = κ[4, :]
    end

    close(f)

    #@info "File $(path) read."
    MARCSOS(T, pe, pg, ρ, κ_c, κ_s, κ_la, κ_lm, λ, ω, Dict(species[i] => abunds[i] for i in eachindex(species)))
end




#= Conversion functions =#

"""
Interolate a MARCS OS from the raw output to a rho-T grid.
It interpolates on a fixed T grid. 
Note that the T grid will not yet be uniform. The output is a 2D Matrix
for all variables except T and rho.
The results from this routine can be appended to each other if the T values
are mutually exclusive. (they should be)
"""
function interpolate_lnRho(m::MARCSOS, new_axis)
    T    = m.T |> unique 
    nT   = T   |> length
    nPg  = floor(Int, length(m.T) / nT)
    point_assignment = falses(length(m.T), nT)

    T   .= T |> sort
    pe   = zeros(eltype(new_axis), nT, length(new_axis))
    pg   = zeros(eltype(new_axis), nT, length(new_axis))
    κ_c  = zeros(eltype(new_axis), nT, length(new_axis), length(m.λ))
    κ_s  = zeros(eltype(new_axis), nT, length(new_axis), length(m.λ))
    κ_la = zeros(eltype(new_axis), nT, length(new_axis), length(m.λ))
    κ_lm = zeros(eltype(new_axis), nT, length(new_axis), length(m.λ))

    #oldV = zeros(eltype(new_axis), nPg) 
    #oldR = zeros(eltype(new_axis), nPg) 

    @inbounds for i in eachindex(m.T)
        j = findfirst(T .≈ m.T[i])  ## Which grid point is this
        point_assignment[i, j] = true
    end

    # wo walk through every fixed T and interpolate all rho to new grid
    # We save the result as a Matrix
    @inbounds for i in 1:nT
        oldR = m.ρ[view(point_assignment, :, i)] .|> log

        oldV = m.pg[view(point_assignment, :, i)] .|> log
        pg[i, :] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)

        oldV .= view(m.pe, view(point_assignment, :, i)) .|> log
        pe[i, :] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)
        
        @inbounds for l in eachindex(m.λ)    
            oldV .= view(m.κ_c, view(point_assignment, :, i), l) .|> log
            κ_c[i, :, l] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)

            oldV .= view(m.κ_s, view(point_assignment, :, i), l) .|> log
            κ_s[i, :, l] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)

            oldV .= view(m.κ_la, view(point_assignment, :, i), l) .|> log
            κ_la[i, :, l] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)

            oldV .= view(m.κ_lm, view(point_assignment, :, i), l) .|> log
            κ_lm[i, :, l] .= linear_interpolation(oldR, oldV, extrapolation_bc=Line()).(new_axis)
        end
    end
    
    MARCSOS(log.(T), pe, pg, new_axis, κ_c, κ_s, κ_la, κ_lm, m.λ, m.ω, m.abundance)
end


"""
Interolate a combined MARCS OS from the interpolated rho grid to the 
uniform lnT grid.
"""
function interpolate_lnT(m::MARCSOS)
    newT = range(minimum(m.T), maximum(m.T), length=length(m.T)) |> collect

    pe   = similar(m.pe) 
    pg   = similar(m.pe)
    κ_c  = similar(m.κ_c) 
    κ_s  = similar(m.κ_c) 
    κ_la = similar(m.κ_c) 
    κ_lm = similar(m.κ_c) 

    oldV = zeros(eltype(newT), length(m.T)) 
    oldT = m.T 

    for i in eachindex(m.ρ)
        oldV     .= view(m.pe, :, i)
        pe[:, i] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)

        oldV     .= view(m.pg, :, i)
        pg[:, i] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)
        
        for l in eachindex(m.λ)    
            oldV         .= view(m.κ_c, :, i, l)
            κ_c[:, i, l] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)

            oldV         .= view(m.κ_s, :, i, l)
            κ_s[:, i, l] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)

            oldV          .= view(m.κ_la, :, i, l)
            κ_la[:, i, l] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)

            oldV          .= view(m.κ_lm, :, i, l)
            κ_lm[:, i, l] .= linear_interpolation(oldT, oldV, extrapolation_bc=Line()).(newT)
        end
    end
    
    MARCSOS(newT, pe, pg, m.ρ, κ_c, κ_s, κ_la, κ_lm, m.λ, m.ω, m.abundance)
end

"""
Square the given MARCSOS part tables. All of the different tables may have different sizes and grids.
We first interpolate every sub-table to a common density grid. Then we pick a common Temp. grid 
between all tables and interpolate 
"""
function square(mos::MARCSOS...)
    m_rho = []
    common_rho = TSO.common_grid(log, mos..., which=:ρ)

    for m in mos
        append!(m_rho, [TSO.interpolate_lnRho(m, common_rho)])
    end

    ## append together all tables and sort them in T
    m_rho = append(m_rho...)

    ## interpolate the table to uniform lnT grid
    interpolate_lnT(m_rho)
end


"""
    complement(mos, eos; lnRoss=:opacity, lnPg=:opacity, lnNe=:opacity)

Create an EoS table from the inerpolated MARCS OS. Specify what quantities should be taken
from where. Electron density will be computed from the electron pressure.

NOTE: This function has not been tested. It is expected to be slow, not optimizes and 
using a lot of RAM. If it takes too long, the individual loops over eos and lambda 
can be drawn together to make it probably much faster. TBD.
"""
function complement(mos, eos; lnEi=:eos, lnRoss=:opacity, lnPg=:opacity, lnNe=:opacity)
    ## The grid is always taken from the MOS
    newlnT   = mos.T
    newlnRho = mos.ρ

    nt = length(newlnT)
    nr = length(newlnRho)
    T  = eltype(newlnT)

    ## Go though additional quantities and take them from where they should be takes
    ### lnEi is bisected from the EoS
    newlnEi = if lnEi == :eos            
        emin,emax,_,_ = TSO.limits(eos)
        E             = zeros(T, nt, nr)
        lims          = [emin, emax]

        for j in eachindex(newlnRho)
            for i in eachindex(newlnT)
                E[i, j] = TSO.bisect(eos, lnRho=newlnRho[j], lnT=newlnT[i], lnEi=lims)
            end
        end

        E
    else
        error("Given lnEi not available. Take from EoS for now.")
    end

    ### lnRoss can be computed or taken. When it should be computed, do it later.
    newlnRoss = if lnRoss == :eos
        ross = zeros(T, nt, nr)
        rho  = similar(newlnT)

        for j in eachindex(newlnRho)
            rho .= newlnRho[j]
            ross[:, j] = TSO.lookup(eos, :lnRoss, rho, newlnT)
        end

        ross
    else
        zeros(T, nt, nr)
    end

    newlnPg = if lnPg == :eos
        y    = zeros(T, nt, nr)
        rho  = similar(newlnT)

        for j in eachindex(newlnRho)
            rho .= newlnRho[j]
            y[:, j] = TSO.lookup(eos, :lnPg, rho, newlnT)
        end

        y
    else
        mos.pg
    end

    newlnNe = if lnNe == :eos
        y    = zeros(T, nt, nr)
        rho  = similar(newlnT)

        for j in eachindex(newlnRho)
            rho .= newlnRho[j]
            y[:, j] = TSO.lookup(eos, :lnNe, rho, newlnT)
        end

        y
    else
        log(exp.(mos.pe) / (TSO.KBoltzmann * exp.(newlnT)))
    end

    # Compute the monochromatic source function (LTE)
    src = zeros(T, nt, nr, length(mos.λ))
    for j in eachindex(mos.λ)
        for i in eachindex(newlnT)
            src[i, :, j] = TSO.Bν(mos.λ[j], exp(newlnT[i]))
        end
    end

    lnT_temp = zeros(T, nt, nr)
    for i in eachindex(newlnRho)
        lnT_temp[:, i] .= newlnT
    end


    ## Combine to EoS + opacity tables
    neweos = TSO.RegularEoSTable(newlnRho, lnT_temp, newlnEi, newlnPg, newlnRoss, newlnNe)

    newopa_c = TSO.RegularOpacityTable(mos.κ_c,                                    newlnRoss, src, mos.λ, false)
    newopa_l = TSO.RegularOpacityTable(mos.κ_lm .+ mos.κ_la,                       newlnRoss, src, mos.λ, false)
    newopa_s = TSO.RegularOpacityTable(mos.κ_s,                                    newlnRoss, src, mos.λ, false)
    #newopa   = TSO.RegularOpacityTable(mos.κ_c .+ mos.κ_lm .+ mos.κ_la .+ mos.κ_s, newlnRoss, src, mos.λ, false)

    ## Now that we have everything we can interpolate to the correct grid (inefficient though)
    neweos, newopa_c, newopa_l, newopa_s = TSO.unify(neweos, [newopa_c, newopa_l, newopa_s])

    ## Recompute the rosseland opacity if wanted
    if lnRoss == :opacity
        TSO.rosseland_opacity!(newlnRoss, neweos, newopa; weights=ω_midpoint(opacities.λ))
        neweos.lnRoss   .= newlnRoss
        newopa_c.κ_ross .= newlnRoss
        newopa_l.κ_ross .= newlnRoss
        newopa_s.κ_ross .= newlnRoss
    end

    newopa   = deepcopy(newopa_c) 
    newopa.κ = newopa_c.κ .+ newopa_l.κ .+ newopa_s.κ

    return neweos, newopa, newopa_c, newopa_l, newopa_s
end

## Utility functions

function common_grid(mos::MARCSOS...; which=:ρ)
    mi = minimum((minimum(getfield(m, which)) for m in mos))
    ma = maximum((maximum(getfield(m, which)) for m in mos))
    l  = maximum((getfield(m, which) |> unique |> length for m in mos))

    range(mi, ma, length=l) |> collect
end

function common_grid(f::Function, mos::MARCSOS...; which=:ρ)
    mi = minimum((minimum(f.(getfield(m, which))) for m in mos))
    ma = maximum((maximum(f.(getfield(m, which))) for m in mos))
    l  = maximum((getfield(m, which) |> unique |> length for m in mos))

    range(mi, ma, length=l) |> collect
end

"""
Append all tables together, sort by increasing T.
"""
function append(mos::TSO.MARCSOS...)
    n_tot = [length(m.T) for m in mos] |> sum
    
    ## Append all the fields
    Tnew    = cat([m.T  for m in mos]..., dims=1)
    penew   = cat([m.pe for m in mos]..., dims=1)
    pgnew   = cat([m.pg for m in mos]..., dims=1)
    ρnew    = mos[1].ρ 
    κ_cnew  = cat([m.κ_c  for m in mos]..., dims=1)
    κ_snew  = cat([m.κ_s  for m in mos]..., dims=1)
    κ_lanew = cat([m.κ_la for m in mos]..., dims=1)
    κ_lmnew = cat([m.κ_lm for m in mos]..., dims=1)

    # Sort them in increasing T order
    mask    = sortperm(Tnew)
    Tnew    = Tnew[mask] 
    penew   = penew[  mask, :]
    pgnew   = pgnew[  mask, :]
    κ_cnew  = κ_cnew[ mask, :, :]
    κ_snew  = κ_snew[ mask, :, :]
    κ_lanew = κ_lanew[mask, :, :]
    κ_lmnew = κ_lmnew[mask, :, :]
    
    MARCSOS(Tnew, penew, pgnew, ρnew, κ_cnew, κ_snew, κ_lanew, κ_lmnew, mos[1].λ, mos[1].ω, mos[1].abundance)
end