# ============================================================================
# EquationOfState -- EoS, opacity and scattering of a table on one grid
# ============================================================================

const EOS_FORMAT = "TSO.EquationOfState"
const EOS_FORMAT_VERSION = 2
const EOS_GRIDS = (:temperature, :energy)
const EOS_PARTS = (:eos, :opacity, :scattering)
const BINNED_EOS_FILENAME = "binned_eos.hdf5"

"""
    EoSData

Thermodynamic quantities of a table, as 2D arrays (energy axis, density axis).
"""
struct EoSData
    data::Dict{Symbol, AbstractArray}
end

EoSData() = EoSData(Dict{Symbol, AbstractArray}())
EoSData(data::AbstractDict) = EoSData(Dict{Symbol, AbstractArray}(Symbol(k) => v for (k, v) in data))

"""
    OpacityData

Opacity-like quantities of a table, as 3D arrays (energy axis, density axis, λ or bin),
together with the wavelength (or bin) axis and its integration weights.
"""
struct OpacityData
    data          ::Dict{Symbol, AbstractArray}
    λ             ::AbstractVector
    weights       ::AbstractVector
    binned        ::Bool
    optical_depth ::Bool
end

_is_binned_axis(λ) = length(λ) > 1 ? all(diff(λ) .== 1) : true

OpacityData(data::AbstractDict, λ; binned=_is_binned_axis(λ), weights=nothing, optical_depth=false) = begin
    d = Dict{Symbol, AbstractArray}(Symbol(k) => v for (k, v) in data)
    w = isnothing(weights) ? (binned ? ones(eltype(λ), length(λ)) : ω_midpoint(to_cm(λ))) : weights
    OpacityData(d, λ, w, binned, optical_depth)
end

OpacityData(κ::AbstractArray{<:Any, 3}, src::Union{Nothing, AbstractArray}, λ::AbstractVector, optical_depth::Bool=false; kwargs...) = begin
    data = Dict{Symbol, AbstractArray}(:κ => κ)
    isnothing(src) || (data[:src] = src)
    OpacityData(data, λ; optical_depth=optical_depth, kwargs...)
end

"""
    EquationOfState

EoS, opacity and scattering of a table on one grid, where the grid is given by
the density axis and either the temperature or the internal energy axis.
Every part (`eos`, `opacity`, `scattering`) may be missing (`nothing`).
All quantities are stored in dictionaries, so any field can be added.

    EquationOfState(path; grid=:temperature, mmap=true, mini=false, parts=(:eos, :opacity, :scattering))

Load a table from a table file or folder. Tables in the old layout (one file per
quantity) are converted on the fly. The opacities are memory-mapped (read-only),
`mmap=false` reads them into memory, which is needed to change them. With
`mini=true` the source functions are not loaded (the Planck function is used
instead) and `parts` selects what is loaded at all.

    tab[:lnPg], tab.lnPg, tab[:κ], tab.λ

Access any quantity of the table, or one of its axes.
"""
struct EquationOfState{D, E, Q<:Union{Nothing, EoSData}, O<:Union{Nothing, OpacityData}, S<:Union{Nothing, OpacityData}}
    density_axis ::D
    energy_axis  ::E
    eos          ::Q
    opacity      ::O
    scattering   ::S
    meta         ::Dict{Symbol, Any}
end

"""
    EquationOfState(lnRho, energy => values; eos=Dict(), opacity=nothing, scattering=nothing, meta=Dict())

Construct a table from its axes. `energy` is `:lnT` or `:lnEi`.
"""
EquationOfState(lnRho::AbstractVector, energy::Pair{Symbol, <:AbstractVector}; eos=EoSData(), opacity=nothing, scattering=nothing, meta=Dict{Symbol, Any}()) = begin
    first(energy) in (:lnT, :lnEi) || error("The energy axis has to be :lnT or :lnEi, got $(first(energy)).")
    e = eos isa AbstractDict ? EoSData(eos) : eos
    EquationOfState(_axis(lnRho, :lnRho, 2), _axis(last(energy), first(energy), 1), e, opacity, scattering, Dict{Symbol, Any}(meta))
end

"""
    EquationOfState(lnRho, lnT, lnEi, lnPg, lnRoss, lnNe; kwargs...)

Construct a table from the standard quantities. The one of `lnT` and `lnEi`
that is a vector is the energy axis.
"""
EquationOfState(lnRho::AbstractVector, lnT::AbstractArray, lnEi::AbstractArray, lnPg::AbstractArray, lnRoss::AbstractArray, lnNe::AbstractArray; kwargs...) = begin
    eos = Dict{Symbol, AbstractArray}(:lnPg => lnPg, :lnRoss => lnRoss, :lnNe => lnNe)
    if ndims(lnEi) == 1
        eos[:lnT] = lnT
        EquationOfState(lnRho, :lnEi => lnEi; eos=eos, kwargs...)
    else
        eos[:lnEi] = lnEi
        EquationOfState(lnRho, :lnT => lnT; eos=eos, kwargs...)
    end
end

"""
    EquationOfState(tab; eos, opacity, scattering, meta)

Copy of `tab` (sharing its data) with some of its parts replaced.
"""
EquationOfState(tab::EquationOfState; eos=tab.eos, opacity=tab.opacity, scattering=tab.scattering, meta=tab.meta) =
    EquationOfState(tab.density_axis, tab.energy_axis, eos, opacity, scattering, Dict{Symbol, Any}(meta))

EquationOfState(opacity::OpacityData; meta=Dict{Symbol, Any}()) =
    EquationOfState(nothing, nothing, nothing, opacity, nothing, Dict{Symbol, Any}(meta))

function _axis(values, name, dimension)
    v = collect(values)
    uni = is_uniform(v)
    EoSAxis(v, dimension, length(v), name, range(minimum(v), maximum(v), length=length(v)), uni)
end

# ============================================================================
# Access
# ============================================================================

Base.getindex(d::EoSData, var::Symbol) = d.data[var]
Base.setindex!(d::EoSData, v, var::Symbol) = (d.data[var] = v)
Base.haskey(d::EoSData, var::Symbol) = haskey(d.data, var)
Base.keys(d::EoSData) = keys(d.data)

Base.getindex(d::OpacityData, var::Symbol) = var == :λ ? d.λ : d.data[var]
Base.setindex!(d::OpacityData, v, var::Symbol) = (d.data[var] = v)
Base.haskey(d::OpacityData, var::Symbol) = var == :λ || haskey(d.data, var)
Base.keys(d::OpacityData) = keys(d.data)
wavelength(d::OpacityData) = d.λ
Base.length(d::OpacityData) = length(d.λ)

has_eos(tab::EquationOfState) = !isnothing(tab.eos)
has_opacity(tab::EquationOfState) = !isnothing(tab.opacity)
has_scattering(tab::EquationOfState) = !isnothing(tab.scattering)
has_axes(tab::EquationOfState) = !isnothing(tab.density_axis)

function Base.getindex(tab::EquationOfState, var::Symbol)
    if has_axes(tab)
        var == :lnRho && return getfield(tab, :density_axis).values
        var == getfield(tab, :energy_axis).name && return getfield(tab, :energy_axis).values
    end
    has_eos(tab) && haskey(getfield(tab, :eos), var) && return getfield(tab, :eos)[var]
    has_opacity(tab) && haskey(getfield(tab, :opacity), var) && return getfield(tab, :opacity)[var]
    error("$(var) not in table, available: $(variables(tab))")
end

function Base.setindex!(tab::EquationOfState, v, var::Symbol)
    if has_opacity(tab) && (haskey(tab.opacity, var) || ndims(v) == 3)
        tab.opacity[var] = v
    elseif has_eos(tab)
        tab.eos[var] = v
    else
        error("No part of the table to store $(var) in.")
    end
end

Base.getproperty(tab::EquationOfState, s::Symbol) = s in fieldnames(EquationOfState) ? getfield(tab, s) : tab[s]
Base.propertynames(tab::EquationOfState) = (fieldnames(EquationOfState)..., variables(tab)...)
Base.haskey(tab::EquationOfState, var::Symbol) = var in variables(tab)

function variables(tab::EquationOfState)
    v = Symbol[]
    has_axes(tab) && push!(v, :lnRho, tab.energy_axis.name)
    has_eos(tab) && append!(v, keys(tab.eos))
    has_opacity(tab) && (push!(v, :λ); append!(v, keys(tab.opacity)))
    unique(v)
end

const EOS_UNITS = Dict{Symbol, String}(
    :lnRho => "ln(g cm-3)", :lnT => "ln(K)", :lnEi => "ln(erg g-1)", :lnPg => "ln(dyn cm-2)",
    :lnNe => "ln(cm-3)", :lnRoss => "ln(cm2 g-1)", :lnK500 => "ln(cm2 g-1)",
    :nan_mask => "boolean", :nan_mask500 => "boolean",
)
const MONOCHROMATIC_UNITS = Dict{Symbol, String}(:κ => "cm2 g-1", :src => "erg s-1 cm-2 sr-1 cm-1", :λ => "Angstrom")
const BINNED_UNITS = Dict{Symbol, String}(:κ => "cm-1", :src => "erg s-1 cm-2 sr-1", :λ => "bin", :ϵ => "dimensionless", :wthin => "dimensionless")

"""
    units(tab)

Units of the quantities of the table. The binned opacity is per cm (it includes
the density) and the binned source function is integrated over the bin, the
monochromatic ones are per gram and per cm of wavelength.
"""
units(tab::EquationOfState) = begin
    binned = has_opacity(tab) ? tab.opacity.binned : get(tab.meta, :kind, "") == "binned"
    u = merge(EOS_UNITS, binned ? BINNED_UNITS : MONOCHROMATIC_UNITS)
    v = variables(tab)
    has_scattering(tab) && append!(v, keys(tab.scattering))
    Dict(k => u[k] for k in v if haskey(u, k))
end

meta(tab::EquationOfState) = tab.meta
meta(tab::EquationOfState, key::Symbol) = tab.meta[key]
table(tab::EquationOfState) = tab
opacity(tab::EquationOfState) = tab.opacity
wavelength(tab::EquationOfState) = wavelength(tab.opacity)
is_binned(tab::EquationOfState) = has_opacity(tab) ? tab.opacity.binned : get(tab.meta, :kind, "") == "binned"
TableInfo(tab::EquationOfState) = get(tab.meta, :table_info, nothing)

"""
    scattering(tab)

The table with the scattering opacity as its opacity.
"""
scattering(tab::EquationOfState) = EquationOfState(tab; opacity=tab.scattering, scattering=nothing)

"""
    merge(eos_tab, opacity_tab, [scattering_tab])

Combine tables of the same grid into one, e.g. after loading files in the old layout separately.
"""
Base.merge(a::EquationOfState, b::EquationOfState, c::Union{Nothing, EquationOfState}=nothing) = EquationOfState(
    has_axes(a) ? a.density_axis : b.density_axis,
    has_axes(a) ? a.energy_axis : b.energy_axis,
    has_eos(a) ? a.eos : b.eos,
    has_opacity(b) ? b.opacity : a.opacity,
    isnothing(c) ? (has_scattering(b) ? b.scattering : a.scattering) : c.opacity,
    merge(a.meta, b.meta, isnothing(c) ? Dict{Symbol, Any}() : c.meta)
)

EnergyAxis(tab::EquationOfState) = tab.energy_axis
DensityAxis(tab::EquationOfState) = tab.density_axis
is_internal_energy(tab::EquationOfState) = tab.energy_axis.name == :lnEi
energy_variable(tab::EquationOfState) = tab.energy_axis.name
dependent_energy_variable(tab::EquationOfState) = is_internal_energy(tab) ? :lnT : :lnEi
grid_name(tab::EquationOfState) = is_internal_energy(tab) ? :energy : :temperature
is_uniform(tab::EquationOfState) = is_uniform(tab.energy_axis) & is_uniform(tab.density_axis)
Base.size(tab::EquationOfState) = (tab.energy_axis.length, tab.density_axis.length)
meshgrid(tab::EquationOfState) = meshgrid(tab.energy_axis.values, tab.density_axis.values)

limits(tab::EquationOfState, which=Colon()) = begin
    e, d = extrema(tab.energy_axis.values), extrema(tab.density_axis.values)
    which == Colon() ? (e..., d...) : (which == 2 ? d : e)
end

Base.show(io::IO, tab::EquationOfState) = begin
    kind = get(tab.meta, :kind, "table")
    parts = [p for p in EOS_PARTS if !isnothing(getfield(tab, p))]
    grid = has_axes(tab) ? "$(grid_name(tab)), size=$(size(tab))" : "no axes"
    nλ = has_opacity(tab) ? ", nλ=$(length(tab.opacity))" : ""
    print(io, "EquationOfState($(kind), $(grid)$(nλ), parts=$(parts))")
end

Base.show(io::IO, ::MIME"text/plain", tab::EquationOfState) = begin
    show(io, tab)
    println(io)
    has_eos(tab) && println(io, "  eos:        ", join(sort(String.(collect(keys(tab.eos)))), ", "))
    has_opacity(tab) && println(io, "  opacity:    ", join(sort(String.(collect(keys(tab.opacity)))), ", "))
    has_scattering(tab) && println(io, "  scattering: ", join(sort(String.(collect(keys(tab.scattering)))), ", "))
    info = TableInfo(tab)
    isnothing(info) || print(io, "  ", show_composition(info))
end

# ============================================================================
# Finding tables
# ============================================================================

"""
    is_eos_file(path)

Check if `path` is a table in the single-file layout (by its attributes, not its name).
"""
is_eos_file(path::AbstractString) = begin
    (isfile(path) && endswith(path, ".hdf5") && HDF5.ishdf5(path)) || return false
    HDF5.h5open(path, "r") do fid
        attrs = HDF5.attrs(fid)
        haskey(attrs, "format") && attrs["format"] == EOS_FORMAT
    end
end

"""
    eos_file(folder)

The single-file table in `folder`, or `nothing` if the folder is in the old layout.
"""
function eos_file(folder::AbstractString)
    isdir(folder) || return nothing
    candidates = filter(is_eos_file, readdir(folder, join=true))
    length(candidates) > 1 && error("More than one table file in $(folder): $(basename.(candidates))")
    isempty(candidates) ? nothing : first(candidates)
end

"""
    availableGrids(path)

Grids (`:temperature`, `:energy`) that are available for the table at `path`.
"""
availableGrids(tab::EquationOfState) = availableGrids(tab.meta[:path])
function availableGrids(path::AbstractString)
    file = isdir(path) ? eos_file(path) : path
    if !isnothing(file) && is_eos_file(file)
        HDF5.h5open(file, "r") do fid
            [g for g in EOS_GRIDS if haskey(fid, String(g))]
        end
    else
        folder = isdir(path) ? path : dirname(path)
        [g for g in EOS_GRIDS if haskey(_legacy_files(folder), g)]
    end
end

# ============================================================================
# Reading
# ============================================================================

function EquationOfState(path::AbstractString; grid=nothing, mmap=true, mini=false, parts=EOS_PARTS)
    path = abspath(path)
    file = isdir(path) ? eos_file(path) : (is_eos_file(path) ? path : nothing)
    skip = mini ? (:src,) : ()

    tab = if !isnothing(file)
        _read_eos_file(file, grid; mmap=mmap, skip=skip, parts=parts)
    elseif isdir(path)
        _read_legacy_folder(path, grid; mmap=mmap, skip=skip, parts=parts)
    elseif isfile(path)
        _read_legacy_file(path; mmap=mmap, skip=skip)
    else
        error("No table found at $(path).")
    end

    tab.meta[:path] = isnothing(file) ? path : file
    tab
end

"""
    switch(tab; grid)

Load the same table on the other grid.
"""
switch(tab::EquationOfState; grid, mmap=true) = EquationOfState(tab.meta[:path]; grid=grid, mmap=mmap)

"""
    reload(path; mmap=true)

Load any table file, in the old or in the new layout, as `EquationOfState`.
"""
reload(path::AbstractString; kwargs...) = EquationOfState(path; kwargs...)

function _read_eos_file(file, grid; mmap=false, skip=(), parts=EOS_PARTS)
    HDF5.h5open(file, "r") do fid
        available = [g for g in EOS_GRIDS if haskey(fid, String(g))]
        isempty(available) && error("$(file) contains no grid.")
        g = isnothing(grid) ? first(available) : Symbol(grid)
        g in available || error("Grid $(g) not in $(file), available: $(available)")

        group = fid[String(g)]
        energy = only(n for n in ("lnT", "lnEi") if haskey(group, n))
        meta = _read_meta(fid)
        binned = get(meta, :kind, "") == "binned"

        has(p) = p in parts && haskey(group, String(p))
        EquationOfState(
            HDF5.read(group["lnRho"]), Symbol(energy) => HDF5.read(group[energy]);
            eos=has(:eos) ? EoSData(_read_group(group["eos"])) : nothing,
            opacity=has(:opacity) ? _read_opacity_group(group["opacity"]; mmap=mmap, binned=binned, skip=skip) : nothing,
            scattering=has(:scattering) ? _read_opacity_group(group["scattering"]; mmap=mmap, binned=binned, skip=skip) : nothing,
            meta=meta
        )
    end
end

_read_array(d, mmap) = (mmap && ndims(d) > 2 && HDF5.ismmappable(d)) ? HDF5.readmmap(d) : HDF5.read(d)

_read_group(group; mmap=false, skip=()) = Dict{Symbol, AbstractArray}(
    Symbol(k) => _read_array(group[k], mmap) for k in keys(group) if !(Symbol(k) in skip)
)

const _DROPPED_FIELDS = (:κ_ross,)

function _read_opacity_group(group; mmap=false, binned=nothing, skip=())
    data = _read_group(group; mmap=mmap, skip=(skip..., _DROPPED_FIELDS...))
    λ = pop!(data, :λ)
    attrs = HDF5.attrs(group)
    optical_depth = haskey(attrs, "optical_depth") ? Bool(attrs["optical_depth"]) : false
    OpacityData(data, λ; binned=isnothing(binned) ? _is_binned_axis(λ) : binned, optical_depth=optical_depth)
end

function _read_meta(fid)
    meta = Dict{Symbol, Any}()
    for (k, v) in pairs(HDF5.attrs(fid))
        k in ("format", "format_version", "table_info") && continue
        meta[Symbol(k)] = v
    end
    info = _read_table_info_h5(fid)
    isnothing(info) || (meta[:table_info] = info)
    meta
end

# ============================================================================
# Writing
# ============================================================================

"""
    save(tab::EquationOfState, path; mode="w")

Save `tab` as the group of its grid in the file `path`. With `mode="r+"` the
grid is added to (or replaced in) an existing file, the other grids are kept.
The meta information is written as attributes of the file.
"""
function save(tab::EquationOfState, path::AbstractString; mode="w")
    has_axes(tab) || error("A table without axes can not be saved, merge it with its EoS first.")
    mode = (mode == "r+" && !isfile(path)) ? "w" : mode
    HDF5.h5open(path, mode) do fid
        _write_meta!(fid, tab.meta)
        g = String(grid_name(tab))
        haskey(fid, g) && HDF5.delete_object(fid, g)
        group = HDF5.create_group(fid, g)
        u = units(tab)
        _write_dataset!(group, :lnRho, tab.density_axis.values, u)
        _write_dataset!(group, tab.energy_axis.name, tab.energy_axis.values, u)
        has_eos(tab) && _write_group!(HDF5.create_group(group, "eos"), tab.eos.data, u)
        has_opacity(tab) && _write_opacity_group!(HDF5.create_group(group, "opacity"), tab.opacity, u)
        has_scattering(tab) && _write_opacity_group!(HDF5.create_group(group, "scattering"), tab.scattering, u)
    end
    path
end

"""
    save(path, tabs::EquationOfState...)

Save several grids of the same table in one file.
"""
save(path::AbstractString, tab::EquationOfState, tabs::EquationOfState...) = begin
    save(tab, path)
    for t in tabs
        save(t, path; mode="r+")
    end
    path
end

function _write_group!(group, data, u)
    for (k, v) in data
        _write_dataset!(group, k, v isa BitArray ? collect(v) : v, u)
    end
end

function _write_opacity_group!(group, opa::OpacityData, u)
    _write_group!(group, opa.data, u)
    _write_dataset!(group, :λ, collect(opa.λ), u)
    HDF5.attrs(group)["optical_depth"] = Int(opa.optical_depth)
end

function _write_meta!(fid, meta)
    attrs = HDF5.attrs(fid)
    for (k, v) in [("format", EOS_FORMAT), ("format_version", EOS_FORMAT_VERSION), ((String(k), v) for (k, v) in meta if !(k in (:path, :table_info)) && !isnothing(v))...]
        haskey(attrs, k) && HDF5.delete_attribute(fid, k)
        attrs[k] = _attribute_value(v)
    end
    isnothing(get(meta, :table_info, nothing)) || _write_table_info_h5!(fid, meta[:table_info])
end

_upgrade_meta!(file) = HDF5.h5open(fid -> _write_meta!(fid, _read_meta(fid)), file, "r+")
_format_version(file) = HDF5.h5open(fid -> Int(get(HDF5.attrs(fid), "format_version", 1)), file, "r")

_write_dataset!(group, name, v, u) = begin
    group[String(name)] = v
    _write_unit!(group[String(name)], name, u)
end

_write_unit!(dataset, name, u) = haskey(u, Symbol(name)) && (HDF5.attrs(dataset)["unit"] = u[Symbol(name)])

_attribute_value(v::Bool) = Int(v)
_attribute_value(v::Symbol) = String(v)
_attribute_value(v) = v

# ============================================================================
# Old layout: one file per quantity and grid
# ============================================================================

const _LEGACY_BINNED = (
    (eos="eos_T.hdf5", eos500="eos_T500.hdf5", opacity="binned_opacities_T.hdf5"),
    (eos="eos.hdf5",   eos500=nothing,         opacity="binned_opacities.hdf5"),
)

_legacy_grid(path) = HDF5.h5open(fid -> ndims(fid["lnEi"]) == 1 ? :energy : :temperature, path, "r")

function _legacy_monochromatic(folder)
    files = readdir(folder)
    pick(prefix) = begin
        f = filter(n -> startswith(n, prefix) && endswith(n, ".hdf5"), files)
        isempty(f) ? nothing : joinpath(folder, first(f))
    end
    eos = pick("combined_eos_")
    isnothing(eos) && return nothing
    (eos=eos, eos500=pick("combined_eos500_"), opacity=pick("combined_opacities_"), scattering=pick("combined_sopacities_"))
end

function _legacy_files(folder)
    found = Dict{Symbol, NamedTuple}()
    mono = _legacy_monochromatic(folder)
    if !isnothing(mono)
        found[:temperature] = mono
        return found
    end
    p(n) = (!isnothing(n) && isfile(joinpath(folder, n))) ? joinpath(folder, n) : nothing
    for names in _LEGACY_BINNED
        eos = p(names.eos)
        isnothing(eos) && continue
        g = _legacy_grid(eos)
        haskey(found, g) && error("$(basename(found[g].eos)) and $(names.eos) in $(folder) are both on the $(g) grid.")
        found[g] = (eos=eos, eos500=p(names.eos500), opacity=p(names.opacity), scattering=nothing)
    end
    found
end

_legacy_kind(folder) = isnothing(_legacy_monochromatic(folder)) ? "binned" : "monochromatic"

function _legacy_meta(folder)
    meta = Dict{Symbol, Any}(:kind => _legacy_kind(folder))
    has_table_info(folder) && (meta[:table_info] = _read_table_info(joinpath(folder, TABLE_INFO_FILENAME)))
    meta
end

_is_legacy_eos(fid) = haskey(fid, "lnEi") && haskey(fid, "lnRho")

function _read_legacy_eos(path)
    data = HDF5.h5open(path, "r") do fid
        Dict{Symbol, AbstractArray}(Symbol(k) => HDF5.read(fid[k]) for k in keys(fid))
    end
    lnRho = pop!(data, :lnRho)
    energy = ndims(data[:lnEi]) == 1 ? :lnEi : :lnT
    EquationOfState(lnRho, energy => pop!(data, energy); eos=data)
end

function _read_legacy_opacity(path; mmap=false, binned=nothing, skip=())
    HDF5.h5open(path, "r") do fid
        data = Dict{Symbol, AbstractArray}()
        for k in keys(fid)
            (k in ("λ", "optical_depth") || Symbol(k) in skip || Symbol(k) in _DROPPED_FIELDS) && continue
            data[Symbol(k)] = _read_array(fid[k], mmap)
        end
        λ = HDF5.read(fid["λ"])
        optical_depth = haskey(fid, "optical_depth") ? Bool(HDF5.read(fid["optical_depth"])) : false
        OpacityData(data, λ; binned=isnothing(binned) ? _is_binned_axis(λ) : binned, optical_depth=optical_depth)
    end
end

function _read_legacy_folder(folder, grid; mmap=false, skip=(), parts=EOS_PARTS)
    files = _legacy_files(folder)
    isempty(files) && error("No table found in $(folder).")
    g = isnothing(grid) ? (haskey(files, :temperature) ? :temperature : first(keys(files))) : Symbol(grid)
    haskey(files, g) || error("Grid $(g) not in $(folder), available: $(collect(keys(files)))")
    f = files[g]
    binned = _legacy_kind(folder) == "binned"

    tab = _read_legacy_eos(f.eos)
    if !isnothing(f.eos500)
        tab500 = _read_legacy_eos(f.eos500)
        tab.eos[:lnK500] = tab500.eos[:lnRoss]
        haskey(tab500.eos, :nan_mask) && (tab.eos[:nan_mask500] = tab500.eos[:nan_mask])
    end

    EquationOfState(tab;
        eos=:eos in parts ? tab.eos : nothing,
        opacity=(:opacity in parts && !isnothing(f.opacity)) ? _read_legacy_opacity(f.opacity; mmap=mmap, binned=binned, skip=skip) : nothing,
        scattering=(:scattering in parts && !isnothing(f.scattering)) ? _read_legacy_opacity(f.scattering; mmap=mmap, binned=binned, skip=skip) : nothing,
        meta=_legacy_meta(folder)
    )
end

function _read_legacy_file(path; mmap=false, skip=())
    is_eos = HDF5.h5open(_is_legacy_eos, path, "r")
    is_eos ? _read_legacy_eos(path) : EquationOfState(_read_legacy_opacity(path; mmap=mmap, skip=skip))
end

# ============================================================================
# Migration
# ============================================================================

"""
    eos_filename(folder)

Name of the single-file table for `folder`: `binned_eos.hdf5` for binned tables,
the folder name for monochromatic tables.
"""
eos_filename(folder) = _legacy_kind(folder) == "binned" ? BINNED_EOS_FILENAME : basename(rstrip(abspath(folder), '/')) * ".hdf5"

"""
    migrateTable(folder; delete=false)

Combine the old files of the table in `folder` into one file. With `delete=true`
the old files are removed after the new file has been compared to them; this
also works later on, for a folder that was migrated before. The meta information
of a table file in an older format version is brought up to date.
The DISPATCH tables and `table_info.toml` are not touched.
"""
function migrateTable(folder; delete=false)
    folder = abspath(folder)
    files = _legacy_files(folder)
    target = eos_file(folder)
    if isnothing(target)
        isempty(files) && error("No table found in $(folder).")
        target = _write_migrated(folder, files)
    else
        _format_version(target) < EOS_FORMAT_VERSION && _upgrade_meta!(target)
        delete && !isempty(files) && !_verify_migration(folder, target, files) && error("$(target) differs from the old files in $(folder).")
    end
    delete && foreach(rm, unique(filter(!isnothing, [p for f in values(files) for p in values(f)])))
    target
end

function _write_migrated(folder, files)
    need = sum(filesize, unique(filter(!isnothing, [p for f in values(files) for p in values(f)])))
    free = diskstat(folder).available
    need > 0.95 * free && error("Not enough space in $(folder): need $(round(need/1e9, digits=1)) GB, $(round(free/1e9, digits=1)) GB free.")

    target = joinpath(folder, eos_filename(folder))
    tmp = target * ".tmp"
    try
        for (i, g) in enumerate(sort(collect(keys(files)), by=g -> findfirst(==(g), EOS_GRIDS)))
            large = _large_fields(files[g])
            tab = _read_legacy_folder(folder, g; skip=union(values(large)...))
            save(tab, tmp; mode=i == 1 ? "w" : "r+")
            u = merge(units(tab), _legacy_kind(folder) == "binned" ? BINNED_UNITS : MONOCHROMATIC_UNITS)
            HDF5.h5open(tmp, "r+") do fid
                for (part, names) in pairs(large), k in names
                    HDF5.h5open(src -> HDF5.copy_object(src[String(k)], fid["$(g)/$(part)"], String(k)), files[g][part], "r")
                    _write_unit!(fid["$(g)/$(part)/$(k)"], k, u)
                end
            end
        end
        _verify_migration(folder, tmp, files) || error("Migrated table differs from the old files in $(folder).")
        mv(tmp, target)
    catch
        rm(tmp, force=true)
        rethrow()
    end
    target
end

function _verify_migration(folder, file, files)
    for (g, f) in files
        large = _large_fields(f)
        skip = union(values(large)...)
        old = _read_legacy_folder(folder, g; skip=skip)
        new = _read_eos_file(file, g; skip=skip)
        Set(variables(old)) == Set(variables(new)) || return false
        all(isequal(old[v], new[v]) for v in variables(old)) || return false
        has_scattering(old) == has_scattering(new) || return false
        if has_scattering(old)
            s_old, s_new = scattering(old), scattering(new)
            all(isequal(s_old[v], s_new[v]) for v in variables(s_old)) || return false
        end
        old.meta == filter(kv -> kv.first != :path, new.meta) || return false
        for (part, names) in pairs(large), k in names
            same = HDF5.h5open(f[part], "r") do a
                HDF5.h5open(b -> _same_dataset(a[String(k)], b["$(g)/$(part)/$(k)"]), file, "r")
            end
            same || return false
        end
    end
    true
end

_large_fields(f) = (opacity=_large_fields(f.opacity), scattering=_large_fields(f.scattering))
_large_fields(::Nothing) = Symbol[]
_large_fields(path::AbstractString) = HDF5.h5open(fid -> [Symbol(k) for k in keys(fid) if ndims(fid[k]) > 2 && !(Symbol(k) in _DROPPED_FIELDS)], path, "r")

function _same_dataset(a, b; block_bytes=2^28)
    (size(a) == size(b) && eltype(a) == eltype(b)) || return false
    n = max(1, block_bytes ÷ (sizeof(eltype(a)) * size(a, 1) * size(a, 2)))
    all(isequal(a[:, :, r], b[:, :, r]) for r in Iterators.partition(1:size(a, 3), n))
end
