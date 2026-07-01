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

    TableInfo(directory::String)

Load a `TableInfo` from the `table_info.toml` file in the given directory.
If the file does not exist, returns `nothing`.

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
    TableInfo(directory::String)

Load a `TableInfo` from the `table_info.toml` file in `directory`.
Returns `nothing` if the file does not exist.
"""
function TableInfo(directory::String)
    path = joinpath(directory, TABLE_INFO_FILENAME)
    isfile(path) || return nothing
    _read_table_info(path)
end

# ============================================================================
# Reading
# ============================================================================

function _read_table_info(path::String)
    data = TOML.parsefile(path)

    feh     = Float64(get(data, "feh",     0.0))
    alpha   = Float64(get(data, "alpha",   0.0))
    vmic    = Float64(get(data, "vmic",    1.0))
    version = String(get(data,  "version", ""))

    # composition is stored as string "Mn_-0.1,Co_0.1" or as TOML sub-table
    raw_comp = get(data, "composition", "")
    composition = if raw_comp isa AbstractString
        parse_composition(raw_comp)
    elseif raw_comp isa Dict
        Dict{String, Float64}(string(k) => Float64(v) for (k, v) in raw_comp)
    else
        Dict{String, Float64}()
    end

    # Everything else goes into extras
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

    data = Dict{String, Any}(
        "feh"         => info.feh,
        "alpha"       => info.alpha,
        "composition" => composition_to_string(info.composition),
        "vmic"        => info.vmic,
        "version"     => info.version,
    )

    # Merge extras (extras cannot overwrite core keys)
    for (k, v) in info.extras
        k in _TABLE_INFO_CORE_KEYS && continue
        isnothing(v) && continue
        data[k] = v
    end

    open(path, "w") do io
        println(io, "# EoS / Opacity Table Info")
        println(io, "# Auto-generated by TSO.jl — edit with care")
        println(io)
        TOML.print(io, data)
    end

    path
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
        isdir(entry) || continue
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

has_table_info(directory::String) = isfile(joinpath(directory, TABLE_INFO_FILENAME))

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
