### Loading paths for wrapper and TS itself

const TSwrapper = Ref("")
const TS        = Ref("")

"""Load the path of the TS codes."""
load(ts_help::Ref{String}, path::String) = begin
    ts_help[] = isdir(abspath(path)) ? abspath(path) : ""
    ts_help[]
end

load_wrapper(path=joinpath(@__FILE__, "../../../../TurboSpectrum-Wrapper") ) = load(TSwrapper, path)
load_TS(path=joinpath(@__FILE__, "../../../../Turbospectrum_NLTE") )         = load(TS, path)


### General Convenience

"""The given path is converted from relative to the wrapper to absolute."""
macro inWrapper(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_wrapper($(relative_path_l)))
end

"""The given path is converted from relative to the TS root dir to absolute."""
macro inTS(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_ts($(relative_path_l)))
end

_in_wrapper(relative_path) = TSwrapper[] == "" ? abspath(joinpath(load_wrapper(), relative_path)) : abspath(joinpath(TSwrapper[], relative_path))
_in_ts(relative_path)      = TS[] == ""        ? abspath(joinpath(load_TS(), relative_path))      : abspath(joinpath(TS[], relative_path))


### Importing Python scripts

function _get_help_py(mod, dir)
	sys   = pyimport("sys")
	dir in sys."path" ? nothing : append!(sys."path",[dir])
    pyimport(String(mod))
end

macro get_help(mod, dir=dirname(@inWrapper("source")))
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end

macro fromWrapper(mod, dir=dirname(@inWrapper("source")))
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end

macro fromTS(mod, dir=dirname(@inTS("source")))
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end


### Writing colunms in Stagger format

"""Write a column to file using the Stagger format as specified in TS"""
function write_as_stagger(name::String; teff, logg, FeH, T, rho, z=zeros(length(T)), id="TSO.jl-column")
    teff = @sprintf "%.1f" teff
    logg = @sprintf "%.5f" logg
    FeH  = @sprintf "%.2f" FeH
    n    = @sprintf "%i" length(z)

    content  = "Stagger  $(id)\n"
    content *= "$(teff) $(logg) $(FeH) $(n)\n"

    for row in eachindex(z)
        content *= @sprintf "%13.5e %13.2f %13.5e\n" z[row] T[row] rho[row]
    end

    open(name, "w") do f
        write(f, content)
    end

    nothing
end

write_as_stagger(name::String, T::Vector, rho::Vector) = write_as_stagger(name, T=T, rho=rho, z=zeros(length(T)), teff=0.0, logg=0.0, FeH=0.0)


### Communication with the TS wrapper
