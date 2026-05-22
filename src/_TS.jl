function write_as_stagger(lnT::Vector, lnRho::Vector; folder=@inWrapper("example/models"), teff=5777.0, logg=4.43, FeH=0.0) 
    sT = length(lnT)
    names_cols = String[]

    nbins = 1#ceil(Int, length(lnT) / 100)
    ibins = TSO.split_similar(lnT, nbins, mask=true)

    # write columns to file
    for column in eachindex(lnRho)
        for i in 1:nbins
            name = joinpath(folder, "TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")
            write_as_stagger(name; teff=teff, logg=logg, FeH=FeH, 
                                    T=exp.(lnT[ibins[i]]), rho=zeros(sT)[ibins[i]] .+ exp(lnRho[column]), id="TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")

            append!(names_cols, ["TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat"])
        end
    end

    open(joinpath(folder, "TSO_list.in"), "w") do f
        for name in names_cols
            write(f, name*"\n")
        end
    end
end

function write_as_stagger(lnT::Matrix, lnRho::Matrix; folder=@inWrapper("example/models"), teff=5777.0, logg=4.43, FeH=0.0) 
    sT = length(lnT)
    names_cols = String[]

    nbins = 1#ceil(Int, size(lnT, 1) / 100)
    ibins = TSO.split_similar(lnT[:, 1], nbins, mask=true)

    # write columns to file
    for column in axes(lnRho, 2)
        for i in 1:nbins
            name = joinpath(folder, "TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")
            write_as_stagger(name; teff=teff, logg=logg, FeH=FeH, 
                                    T=exp.(lnT[ibins[i], column]), rho=exp.(lnRho[ibins[i], column]), id="TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat")

            append!(names_cols, ["TSOeos_$(column)-$(first(ibins[i])):$(last(ibins[i])).dat"])
        end
    end

    open(joinpath(folder, "TSO_list.in"), "w") do f
        for name in names_cols
            write(f, name*"\n")
        end
    end
end

"""
Read EoS columns from TS. Order: lnT, lnPe, lnRho, lnPg, lnEi, lnOp
"""
function read_eos_column(path, get_header=false)
    # Read the EoS column
    f = FortranFile(path, convert="big-endian")

    # Fist line with the size of the arrays
    # Order: lnT, lnPe, lnRho, lnPg, lnEi, lnOp
    ntau, xls = read(f, Int32, Float32)
    data = zeros(Float32, ntau, 6)

    for i in 1:ntau
        data[i, :] .= read(f, Float32, Float32, Float32, Float32, Float32, Float32)
    end

    close(f)

    return if get_header
        ntau, xls, data
    else
        data
    end
end

"""
Read Opacity columns from TS.
"""
function read_opacity_column(path)
    f = FortranFile(path, convert="big-endian")

    lheader = 0
    while true
        s = try
            Base.convert(String, read(f, FString{2}))
        catch
            " "
        end
        
        (s[end] !== '*') && break
        lheader += 1
    end

    # Skip the header
    FortranFiles.rewind(f)
    for i in 1:lheader
        read(f, FString{2})
    end

    abund = read(f, (Float32, 92))
    ntau, maxlam, λ_start, λ_end, Δλ = read(f, Int32, Int32, Float64, Float64, Float64)
    xls     = read(f, Float32)
    τ, vmic = read(f, (Float32, ntau), (Float32, ntau))
    κ_ross  = read(f, (Float32, ntau))
    κ_cont  = read(f, (Float32, ntau, maxlam))
    κ_line  = read(f, (Float32, ntau, maxlam))
    κ_scat  = read(f, (Float32, ntau, maxlam))
    close(f)

    opacity = zeros(Float32, ntau, maxlam, 3)
    opacity[:, :, 1] .= κ_cont
    opacity[:, :, 2] .= κ_line
    opacity[:, :, 3] .= κ_scat
    λ = Float32[λ_start:Δλ:λ_end...]

    (opacity, κ_ross, λ, abund)
end

matching_opacity(path) = begin
    new_p = path[first(findfirst("TSOeos", path)):findlast('.', path)-1]
    new_p = new_p[1:findlast('.', new_p)] *"dat.multi"
end

function _read_tables_eos_op(list_of_eos_tables, opacity_folder=@inTS(""); get_individual=false, kwargs...)
    # Get the dimensions
    nrows_per_table, ntables, T = get_dimensions(list_of_eos_tables)

    opacity_path = joinpath(opacity_folder, matching_opacity(list_of_eos_tables[1]))
  
    eaxis = false
    
    _,_,λ,abund = read_opacity_column(opacity_path)

    # Read the data
    lnρ    = zeros(T, ntables)
    lnT2   = zeros(T, nrows_per_table, ntables)
    lnPg   = zeros(T, nrows_per_table, ntables)
    lnEi2  = zeros(T, nrows_per_table, ntables)
    lnNe   = zeros(T, nrows_per_table, ntables)
    ross   = zeros(T, nrows_per_table, ntables)
    lnκ500 = zeros(T, nrows_per_table, ntables)
    κ      = zeros(T, nrows_per_table, ntables, length(λ))
    κc     = zeros(T, nrows_per_table, ntables, length(λ))
    κl     = zeros(T, nrows_per_table, ntables, length(λ))
    κs     = zeros(T, nrows_per_table, ntables, length(λ))
    Sν     = similar(κ)

    for (i,table) in enumerate(list_of_eos_tables)
        column, e_r = get_TSO_index(table)
        eos = read_eos_column(table)                                                   # The columns are ln(T), ln(pe), ln(ρ), ln(Pg), ln(E), ln(κ500)
        lnρ[column]        = eos[1, 3]
        lnT2[e_r, column]   .= eos[:, 1]
        lnPg[e_r, column]   .= eos[:, 4]
        lnEi2[e_r, column]  .= eos[:, 5]
        lnNe[e_r, column]   .= log.(exp.(eos[:, 2]) ./ (KBoltzmann .* exp.(eos[:, 1])))  # p = nkT
        lnκ500[e_r, column] .= eos[:, 6]                                                 # This is not the rosseland opacity! This is the 500nm for now

        opacity, κ_ross, _, _ = read_opacity_column(joinpath(opacity_folder, matching_opacity(table)))
        for j in axes(opacity, 1)
            opc = view(opacity, j, :, 1)
            opl = view(opacity, j, :, 2)
            ops = view(opacity, j, :, 3)

            opc[isnan.(opc) .| (opc .< 0.0)] .= 0.0
            opl[isnan.(opl) .| (opl .< 0.0)] .= 0.0
            ops[isnan.(ops) .| (ops .< 0.0)] .= 0.0
        end

        κ[e_r, column, :]  .= opacity[:, :, 1] .+ opacity[:, :, 2] .+ opacity[:, :, 3]
        κc[e_r, column, :] .= opacity[:, :, 1] 
        κl[e_r, column, :] .=                     opacity[:, :, 2] 
        κs[e_r, column, :] .=                                         opacity[:, :, 3]
        ross[e_r, column]  .= κ_ross
        Sν[e_r, column, :] .= Base.convert.(T, Bν(λ, exp.(lnT2[e_r, column]))) # Planck function (i.e. source function)
    end

    # check which of T/E is the axis and which is lnDependent
    for i in axes(lnT2, 2)
        if lnT2[:, i] != lnT2[:, 1]
            @info "Non-equal T columns detected. Assuming energy axis."
            eaxis = true
            break
        end
    end

    #=mask    = sortperm(lnρ)
    lnρ    .= lnρ[mask]
    lnPg   .= lnPg[:, mask]
    lnEi   .= lnEi[:, mask]
    lnNe   .= lnNe[:, mask]
    ross   .= ross[:, mask]
    κ      .= κ[:, mask, :] 
    lnκ500 .= lnκ500[:, mask]=#

    opacity_tables = if get_individual
        (RegularOpacityTable(κ, ross, Sν, λ, false), RegularOpacityTable(κc, ross, Sν, λ, false), RegularOpacityTable(κl, ross, Sν, λ, false), RegularOpacityTable(κs, ross, Sν, λ, false))
    else
        RegularOpacityTable(κ, ross, Sν, λ, false)
    end

    if eaxis
        (RegularEoSTable(lnρ, lnT2, lnEi2, lnPg, log.(ross), lnNe), opacity_tables)
    else
        (RegularEoSTable(lnρ, lnT2[:, 1], lnEi2, lnPg, log.(ross), lnNe), opacity_tables)
    end
end

function _read_tables_eos(list_of_eos_tables; kwargs...)
    # Get the dimensions
    nrows_per_table, ntables, T = get_dimensions(list_of_eos_tables)

    eaxis = false
    
    # Read the data
    lnρ    = zeros(T, ntables)
    lnT2   = zeros(T, nrows_per_table, ntables)
    lnPg   = zeros(T, nrows_per_table, ntables)
    lnEi2  = zeros(T, nrows_per_table, ntables)
    lnNe   = zeros(T, nrows_per_table, ntables)
    lnκ500 = zeros(T, nrows_per_table, ntables) 

    for (i,table) in enumerate(list_of_eos_tables)
        column, e_r = get_TSO_index(table)
        eos = read_eos_column(table)                                                   # The columns are ln(T), ln(pe), ln(ρ), ln(Pg), ln(E), ln(κ500)
        lnρ[column]          = eos[1, 3]
        lnT2[e_r, column]   .= eos[:, 1]
        lnPg[e_r, column]   .= eos[:, 4]
        lnEi2[e_r, column]  .= eos[:, 5]
        lnNe[e_r, column]   .= log.(exp.(eos[:, 2]) ./ (KBoltzmann .* exp.(eos[:, 1]))) # p = nkT
        lnκ500[e_r, column] .= eos[:, 6]                                                # This is not the rosseland opacity! This is the 500nm for now
    end

    # check which of T/E is the axis and which is lnDependent
    for i in axes(lnT2, 2)
        if lnT2[:, i] != lnT2[:, 1]
            @info "Non-equal T columns detected. Assuming energy axis."
            eaxis = true
            break
        end
    end

    #=mask    = sortperm(lnρ)
    lnρ    .= lnρ[mask]
    lnPg   .= lnPg[:, mask]
    lnEi   .= lnEi[:, mask]
    lnNe   .= lnNe[:, mask]
    ross   .= ross[:, mask]
    κ      .= κ[:, mask, :] 
    lnκ500 .= lnκ500[:, mask]=#

    if eaxis
        RegularEoSTable(lnρ, lnT2, lnEi2, lnPg, lnκ500, lnNe)
    else
        RegularEoSTable(lnρ, lnT2[:, 1], lnEi2, lnPg, lnκ500, lnNe)
    end
end

function get_dimensions(list_of_eos_tables)
    e = []
    d = []
    for i in eachindex(list_of_eos_tables)
        column, e_r = get_TSO_index(list_of_eos_tables[i])
        if column == 1
            append!(e, last(e_r))
        end

        append!(d, column)
    end

    eos = read_eos_column(list_of_eos_tables[1])
    #@show maximum(e), length(unique(d)), eltype(eos)
    return maximum(e), length(unique(d)), eltype(eos)
end

function read_tables(list_of_eos_tables::AbstractVector; kwargs...)
    opacity_path = @inTS "TSOeos_$(get_TSO_index(list_of_eos_tables[1])).multi"
    return if ispath(opacity_path)
        _read_tables_eos_op(list_of_eos_tables; kwargs...)
    else
        _read_tables_eos(list_of_eos_tables; kwargs...)
    end
end

read_tables(::Type{<:EoSTable}, ::Type{<:OpacityTable}, list_of_eos_tables::AbstractVector, opacity_folder=@inTS(""); kwargs...) = _read_tables_eos_op(list_of_eos_tables, opacity_folder; kwargs...)
read_tables(::Type{<:EoSTable}, list_of_eos_tables::AbstractVector; kwargs...) = _read_tables_eos(list_of_eos_tables; kwargs...)

read_tables(t::Type{<:EoSTable}, folder::String=".", args...; kwargs...) = read_tables(t, glob("_TSOeos_*_TSO.eos", folder), args...; kwargs...)
read_tables(t::Type{<:EoSTable}, t2::Type{<:OpacityTable}, folder::String=".", args...; kwargs...) = read_tables(t, t2, glob("_TSOeos_*_TSO.eos", folder), args...; kwargs...)

# User interface function
load(e::Type{<:EoSTable}, o::Type{<:OpacityTable}, args...; kwargs...) = read_tables(e, o, args...; kwargs...)
load(e::Type{<:EoSTable}, args...; kwargs...) = read_tables(e,    args...; kwargs...)
load(o::Type{<:OpacityTable}, e::Type{<:EoSTable}, args...; kwargs...) = begin
    e, o = read_tables(e, o, args...; kwargs...)
    o, e
end

get_TSO_index(name) = begin
    mask  = first(findfirst("TSOeos", name))
    maske = first(findlast(".", name)) -1
    rest  = name[mask:maske]
    rest  = rest[1:first(findlast(".", rest))-1]
    ident = split(split(rest, "_")[2], "-")
    d_col = parse(Int, ident[1])
    split(ident[2], ":")
    e_r   = range(parse(Int, split(ident[2], ":")[1]),parse(Int, split(ident[2], ":")[2]), step=1)

    d_col, e_r
end