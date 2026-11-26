import TOML

const TABLE_INFO_FILENAME = "table_info.toml"
const _TABLE_INFO_CORE_KEYS = Set(["feh", "alpha", "composition", "vmic", "version"])

"""
    TableInfo

Metadata for an EoS / opacity table, stored as a TOML file in the table directory.
Replaces the convention of encoding all parameters in the directory name.

# Core fields
- `feh::Float64`                    — Metallicity [Fe/H] (solar-scaled)
- `alpha::Float64`                  — Alpha-element enhancement [α/Fe]
- `composition::Dict{String,Float64}` — Individual element abundances [X/Fe], e.g. `Dict("Mn" => -0.1, "Co" => 0.1)`
- `vmic::Float64`                   — Micro-turbulence velocity [km/s]
- `version::String`                 — Version identifier string (e.g. "v3.6")
- `extras::Dict{String,Any}`        — Any additional key-value pairs found in the info file

# Constructors

    TableInfo(; feh, alpha, composition, vmic, version, kwargs...)

Create a `TableInfo` from keyword arguments. `composition` can be a `Dict`,
or a string in the `"Mn_-0.1,Co_0.1"` format (which is auto-parsed).
Extra keyword arguments are stored in `extras`.

    TableInfo(path::String)

Load the `TableInfo` of a table file or table folder (from the meta information of the
table file, or from `table_info.toml` for tables in the old layout), `nothing` if there is none.

# Examples

```julia
# Create from keyword arguments
info = TableInfo(feh=-1.0, alpha=0.4, composition=Dict("Mn" => -0.1), vmic=1.0, version="v3.6")

# Create from composition string (as passed on the CLI)
info = TableInfo(feh=-1.0, alpha=0.4, composition="Mn_-0.1,Co_0.1", vmic=1.0, version="v3.6")

# Save and reload
save(info, "/path/to/table_dir")
info2 = TableInfo("/path/to/table_dir")

# Access extras
info = TableInfo(feh=0.0, alpha=0.0, vmic=1.0, version="v3.6", scattering=true)
info.extras["scattering"]  # true
```
"""
struct TableInfo
    feh         ::Float64
    alpha       ::Float64
    composition ::Dict{String, Float64}
    vmic        ::Float64
    version     ::String
    extras      ::Dict{String, Any}
end

# ============================================================================
# Composition string <-> Dict conversion
# ============================================================================

function parse_composition(comp_str::AbstractString)
    s = strip(comp_str)
    isempty(s) && return Dict{String, Float64}()
    elements_abundances = split(s, ",", keepempty=false)
    Dict{String, Float64}(
        String(first(split(a, "_", keepempty=false))) => parse(Float64, last(split(a, "_", keepempty=false)))
        for a in elements_abundances
    )
end
parse_composition(comp::AbstractDict) = Dict{String, Float64}(string(k) => Float64(v) for (k, v) in comp)
parse_composition(::Nothing) = Dict{String, Float64}()

function composition_to_string(comp::Dict{String, Float64})
    isempty(comp) && return ""
    join(["$(k)_$(v)" for (k, v) in sort(collect(comp), by=first)], ",")
end

# ============================================================================
# Constructors
# ============================================================================

"""
    TableInfo(; feh=0.0, alpha=0.0, composition=Dict{String,Float64}(), vmic=1.0, version="", kwargs...)

Keyword constructor. `composition` accepts either a `Dict{String,Float64}` or
a composition string (e.g. `"Mn_-0.1,Co_0.1"`) that will be auto-parsed.
Any keyword arguments not in the core set are stored in `extras`.
"""
function TableInfo(; feh::Real=0.0,
                     alpha::Real=0.0,
                     composition=Dict{String, Float64}(),
                     vmic::Real=1.0,
                     version::AbstractString="",
                     kwargs...)
    comp = if composition isa AbstractString
        parse_composition(composition)
    elseif composition isa Dict{Symbol, Float64}
        Dict{String, Float64}(string(k) => v for (k, v) in composition)
    else
        Dict{String, Float64}(string(k) => Float64(v) for (k, v) in composition)
    end
    extras = Dict{String, Any}(string(k) => v for (k, v) in kwargs)
    TableInfo(Float64(feh), Float64(alpha), comp,
              Float64(vmic), String(version), extras)
end

"""
    TableInfo(path::String)

Load the `TableInfo` of a table file, or of the table in the folder `path`, from the
meta information of the table file or, for tables in the old layout, from `table_info.toml`.
Returns `nothing` if neither exists.
"""
function TableInfo(directory::String)
    file = isdir(directory) ? eos_file(directory) : (is_eos_file(directory) ? directory : nothing)
    if !isnothing(file)
        info = HDF5.h5open(_read_table_info_h5, file, "r")
        isnothing(info) || return info
    end
    isdir(directory) || return nothing
    path = joinpath(directory, TABLE_INFO_FILENAME)
    isfile(path) || return nothing
    _read_table_info(path)
end

# ============================================================================
# Reading
# ============================================================================

_read_table_info(path::String) = _table_info(TOML.parsefile(path))
_table_info_from_toml(toml::AbstractString) = _table_info(TOML.parse(toml))
_table_info_to_toml(info::TableInfo) = sprint(TOML.print, _table_info_dict(info))

function _read_table_info_h5(fid)
    haskey(fid, "table_info") && return _table_info(_read_attribute_tree(fid["table_info"]))
    attrs = HDF5.attrs(fid)
    haskey(attrs, "table_info") ? _table_info_from_toml(attrs["table_info"]) : nothing
end

_read_attribute_tree(group) = merge(
    Dict{String, Any}(k => v for (k, v) in pairs(HDF5.attrs(group))),
    Dict{String, Any}(k => _read_attribute_tree(group[k]) for k in keys(group))
)

function _write_table_info_h5!(fid, info::TableInfo)
    haskey(fid, "table_info") && HDF5.delete_object(fid, "table_info")
    haskey(HDF5.attrs(fid), "table_info") && HDF5.delete_attribute(fid, "table_info")
    data = _table_info_dict(info)
    data["composition"] = info.composition
    _write_attribute_tree!(HDF5.create_group(fid, "table_info"), data)
end

function _write_attribute_tree!(group, data)
    for (k, v) in data
        if v isa AbstractDict
            _write_attribute_tree!(HDF5.create_group(group, String(k)), v)
        else
            HDF5.attrs(group)[String(k)] = _attribute_tree_value(v)
        end
    end
end

_attribute_tree_value(v::AbstractVector) = isempty(v) ? Float64[] : collect(promote(v...))
_attribute_tree_value(v) = v

function _table_info(data::AbstractDict)

    feh     = Float64(get(data, "feh",     0.0))
    alpha   = Float64(get(data, "alpha",   0.0))
    vmic    = Float64(get(data, "vmic",    1.0))
    version = String(get(data,  "version", ""))

    raw_comp = get(data, "composition", "")
    composition = if raw_comp isa AbstractString
        parse_composition(raw_comp)
    elseif raw_comp isa Dict
        Dict{String, Float64}(string(k) => Float64(v) for (k, v) in raw_comp)
    else
        Dict{String, Float64}()
    end

    extras = Dict{String, Any}()
    for (k, v) in data
        k in _TABLE_INFO_CORE_KEYS && continue
        extras[k] = v
    end

    TableInfo(feh, alpha, composition, vmic, version, extras)
end

# ============================================================================
# Writing / Saving
# ============================================================================

"""
    save(info::TableInfo, directory::String)

Write the `TableInfo` as a TOML file (`table_info.toml`) in the given directory.
The directory is created if it does not already exist.
The `composition` dictionary is serialised as a string (e.g. `"Co_0.1,Mn_-0.1"`).
Returns the full path to the written file.
"""
function save(info::TableInfo, directory::String)
    isdir(directory) || mkpath(directory)
    path = joinpath(directory, TABLE_INFO_FILENAME)
    data = _table_info_dict(info)

    open(path, "w") do io
        println(io, "# EoS / Opacity Table Info")
        println(io, "# Auto-generated by TSO.jl — edit with care")
        println(io)
        TOML.print(io, data)
    end

    path
end

function _table_info_dict(info::TableInfo)
    data = Dict{String, Any}(
        "feh"         => info.feh,
        "alpha"       => info.alpha,
        "composition" => composition_to_string(info.composition),
        "vmic"        => info.vmic,
        "version"     => info.version,
    )

    for (k, v) in info.extras
        k in _TABLE_INFO_CORE_KEYS && continue
        isnothing(v) && continue
        data[k] = v
    end

    data
end

# ============================================================================
# Matching / Comparison
# ============================================================================

function matches(info::TableInfo;
                 feh::Real=0.0,
                 alpha::Real=0.0,
                 composition=Dict{String, Float64}(),
                 vmic::Real=1.0,
                 version::String="")
    comp = if composition isa AbstractString
        parse_composition(composition)
    elseif composition isa Dict{Symbol, Float64}
        Dict{String, Float64}(string(k) => v for (k, v) in composition)
    else
        Dict{String, Float64}(string(k) => Float64(v) for (k, v) in composition)
    end

    v_match = isempty(version) ? true : (info.version == version)

    v_match &&
    isapprox(info.feh, Float64(feh), atol=1e-6) &&
    isapprox(info.alpha, Float64(alpha), atol=1e-6) &&
    isapprox(info.vmic, Float64(vmic), atol=1e-6) &&
    _compositions_match(info.composition, comp)
end

function _compositions_match(a::Dict{String, Float64}, b::Dict{String, Float64})
    keys(a) == keys(b) || return false
    all(isapprox(a[k], b[k], atol=1e-6) for k in keys(a))
end

function find_table_info(tables_dir::String;
                         feh::Real=0.0,
                         alpha::Real=0.0,
                         composition=Dict{String, Float64}(),
                         vmic::Real=1.0,
                         version::String="")
    isdir(tables_dir) || return nothing

    matches_found = String[]
    for entry in readdir(tables_dir, join=true)
        (isdir(entry) || is_eos_file(entry)) || continue
        info = TableInfo(entry)
        isnothing(info) && continue
        if matches(info; feh=feh, alpha=alpha, composition=composition, vmic=vmic, version=version)
            push!(matches_found, entry)
        end
    end

    if length(matches_found) > 1
        @warn "Multiple tables found matching the same configuration and version:\n" * join(matches_found, "\n")
        error("Duplicate tables detected! The core cannot decide which one to use. Please remove duplicates and try again.")
    elseif length(matches_found) == 1
        return matches_found[1]
    else
        return nothing
    end
end

# ============================================================================
# Utility
# ============================================================================

has_table_info(path::String) = (isdir(path) && isfile(joinpath(path, TABLE_INFO_FILENAME))) || !isnothing(TableInfo(path))

function Base.show(io::IO, info::TableInfo)
    comp_str = composition_to_string(info.composition)
    print(io, "TableInfo(")
    print(io, "feh=$(info.feh), ")
    print(io, "alpha=$(info.alpha), ")
    print(io, "composition=\"$(comp_str)\", ")
    print(io, "vmic=$(info.vmic), ")
    print(io, "version=\"$(info.version)\"")
    if !isempty(info.extras)
        for (k, v) in info.extras
            print(io, ", $(k)=$(repr(v))")
        end
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", info::TableInfo)
    println(io, "TableInfo:")
    println(io, "  [Fe/H]       = $(info.feh)")
    println(io, "  [α/Fe]       = $(info.alpha)")
    comp_str = composition_to_string(info.composition)
    println(io, "  composition  = $(isempty(comp_str) ? "(solar)" : comp_str)")
    println(io, "  vmic [km/s]  = $(info.vmic)")
    println(io, "  version      = $(info.version)")
    if !isempty(info.extras)
        println(io, "  extras:")
        for (k, v) in sort(collect(info.extras), by=first)
            println(io, "    $(k) = $(repr(v))")
        end
    end
end

Base.:(==)(a::TableInfo, b::TableInfo) = (
    a.feh == b.feh &&
    a.alpha == b.alpha &&
    a.composition == b.composition &&
    a.vmic == b.vmic &&
    a.version == b.version &&
    a.extras == b.extras
)

function show_composition(feh, alpha, comp_dict)
    cs = isempty(comp_dict) ? "" : ", with: " * join(["[$(k)/Fe]=$(v)" for (k, v) in comp_dict], ",")
    return "[Fe/H]=$(feh), [α/Fe]=$(alpha)$(cs)"
end

show_composition(info::TableInfo) = show_composition(info.feh, info.alpha, info.composition)

"""
    composition_dict(info)

The composition of a table as `Dict` of `Symbol("[Fe/H]")`, `Symbol("[alpha/Fe]")` and
`Symbol("[X/Fe]")` to their value, e.g. for the parameters of a `MUST.Box`.
"""
composition_dict(feh, alpha, composition=Dict()) = Dict{Symbol, Float64}(
    Symbol("[Fe/H]") => feh, Symbol("[alpha/Fe]") => alpha,
    (Symbol("[$(k)/Fe]") => v for (k, v) in parse_composition(composition))...
)
composition_dict(info::TableInfo) = composition_dict(info.feh, info.alpha, info.composition)
composition_dict(::Nothing) = Dict{Symbol, Float64}()
