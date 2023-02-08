## General Utility functions
"""
Split array arr in nsplits roughly equal junks.
If mask=true the indices will be returned.
"""
function split_similar(arr, nsplits; mask=false)
    nsplits = length(arr) < nsplits ? length(arr) : nsplits
    splits  = div(length(arr), nsplits)
    nrest   = length(arr) % nsplits
    split_sizes = Int[splits for _ in 1:nsplits]
    
    for i in 1:nrest
        split_sizes[i] = split_sizes[i] + 1
    end

    split_masks = []
    last_idx    = 0
    for i in 1:nsplits
        append!(split_masks, [[ (last_idx+1:last_idx+split_sizes[i])... ]])
        last_idx = last_idx+split_sizes[i]
    end

    if mask
        return split_masks
    else
        return [arr[mask] for mask in split_masks]
    end
end

function meshgrid(ax...)
    grids = []
    space = Iterators.product(ax...)
    for i in eachindex(ax)
        append!(grids, [getindex.(space, i)])
    end

    grids
end

"""
create random number between a and b
"""
randrange(a,b,args...) = begin
    xmin,xmax = min(a,b),max(a,b)
    rand(args...) .* (xmax-xmin) .+ xmin
end





## Importing of python modules

function _get_help_py(mod ,dir=dirname(@__FILE__))
	sys   = pyimport("sys")
	dir in sys."path" ? nothing : append!(sys."path",[dir])
    pyimport(String(mod))
end

macro pythonHelp(mod)
    mod_e = esc(mod)
    mod_s = :($mod)
    path  = dirname(@__FILE__)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path))
end

macro pythonHelp(mod, dir)
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end




## Saving as HDF5
function save(s::T, path) where {T<:AbstractTable}
    fid = HDF5.h5open(path, "w")

    for fname in fieldnames(typeof(s))
        add_to_hdf5!(fid, fname, getfield(s, fname))
    end

    close(fid)

    path
end

function reload(s::Type{S}, path::String; mmap=false) where {S}
    fid   = HDF5.h5open(path, "r")
    fvals = Any[]

    for (fname, ftype) in zip(fieldnames(s), fieldtypes(s))
        append!(fvals, [get_from_hdf5(ftype, fid, fname; mmap=mmap)])
    end

    close(fid)

    s(fvals...)
end

add_to_hdf5!(fid, fname, val)       = fid["$(fname)"] = val
add_to_hdf5!(fid, fname, val::Bool) = fid["$(fname)"] = Int(val)

get_from_hdf5(::Type{<:Any}, fid, fname; mmap=false) = mmap ? HDF5.readmmap(fid["$(fname)"]) : HDF5.read(fid["$(fname)"])
get_from_hdf5(::Type{Bool},  fid, fname; mmap=false) = Bool(HDF5.read(fid["$(fname)"]))
    
ΔΛ(lo,hi,R)  = (hi+lo)/2 /R
N_Λ(lo,hi,R) = (hi-lo) / ΔΛ(lo,hi,R)





## Constants (In agreement with Tabgen)
const EiExtra    = 5.0                         # extra eV per H atom to add
const KBoltzmann = 1.380658E-16                # Boltzman's cst. [erg/K]
const CLight     = 2.99792458E+10              # Speed of light [cm/s]
const HPlanck    = 6.6260755E-27               # Planck's constant [erg s]
const ev_to_erg  = 1.60218e-12                 # conversion
const HIonPot    = 13.595                      # hydrogen ioniz. potential [eV]
const twohc2     = 2.0e0 *HPlanck*CLight^2
const hc_k       = HPlanck*CLight/KBoltzmann
const aa_to_cm   = 1.0e-8
const σ_S        = 5.6704e-5
const atomic_number = Symbol.(strip(i, ' ') for i in [
        "H ","He","Li","Be","B ","C ","N ","O ","F ",         # |  1 -  9
        "Ne","Na","Mg","Al","Si","P ","S ","Cl","Ar",         # | 10 - 18
        "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co",         # | 19 - 27
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",         # | 28 - 36
        "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh",         # | 37 - 45
        "Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe",         # | 46 - 54
        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu",         # | 55 - 63
        "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",         # | 64 - 72
        "Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl",         # | 73 - 81
        "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",         # | 82 - 90
        "Pa","U "]) 
const atomic_mass = Dict(   "H"  => ("Hydrogen"	    ,1.00797),
                            "He" =>	("Helium"	    ,4.00260),
                            "Li" =>	("Lithium"	    ,6.941),
                            "Be" =>	("Beryllium"	,9.01218),
                            "B"  => ("Boron"	    ,10.81),
                            "C"  => ("Carbon"	    ,12.011),
                            "N"  => ("Nitrogen"	    ,14.0067),
                            "O"  => ("Oxygen"	    ,15.9994),
                            "F"  => ("Fluorine"	    ,18.998403),
                            "Ne" =>	("Neon"	        ,20.179),
                            "Na" =>	("Sodium"	    ,22.98977),
                            "Mg" =>	("Magnesium"	,24.305),
                            "Al" =>	("Aluminum"	    ,26.98154),
                            "Si" =>	("Silicon"	    ,28.0855),
                            "P"  => ("Phosphorus"	,30.97376),
                            "S"  => ("Sulfur"	    ,32.06),
                            "Cl" =>	("Chlorine"	    ,35.453),
                            "K"  => ("Potassium"	,39.0983),
                            "Ar" =>	("Argon"	    ,39.948),
                            "Ca" =>	("Calcium"	    ,40.08),
                            "Sc" =>	("Scandium"	    ,44.9559),
                            "Ti" =>	("Titanium"	    ,47.90),
                            "V"  => ("Vanadium"	    ,50.9415),
                            "Cr" =>	("Chromium"	    ,51.996),
                            "Mn" =>	("Manganese"	,54.9380),
                            "Fe" =>	("Iron"	        ,55.847),
                            "Ni" =>	("Nickel"	    ,58.70),
                            "Co" =>	("Cobalt"	    ,58.9332),
                            "Cu" =>	("Copper"	    ,63.546),
                            "Zn" =>	("Zinc"	        ,65.38),
                            "Ga" =>	("Gallium"	    ,69.72),
                            "Ge" =>	("Germanium"	,72.59),
                            "As" =>	("Arsenic"	    ,74.9216),
                            "Se" =>	("Selenium"	    ,78.96),
                            "Br" =>	("Bromine"	    ,79.904),
                            "Kr" =>	("Krypton"	    ,83.80),
                            "Rb" =>	("Rubidium"	    ,85.4678),
                            "Sr" =>	("Strontium"	,87.62),
                            "Y"  => ( "Yttrium"	    ,88.9059),
                            "Zr" =>	("Zirconium"	,91.22),
                            "Nb" =>	("Niobium"	    ,92.9064),
                            "Mo" =>	("Molybdenum"	,95.94),
                            "Tc" =>	("Technetium"	,98),
                            "Ru" =>	("Ruthenium"	,101.07),
                            "Rh" =>	("Rhodium"	    ,102.9055),
                            "Pd" =>	("Palladium"	,106.4),
                            "Ag" =>	("Silver"	    ,107.868),
                            "Cd" =>	("Cadmium"	    ,112.41),
                            "In" =>	("Indium"	    ,114.82),
                            "Sn" =>	("Tin"	        ,118.69),
                            "Sb" =>	("Antimony"	    ,121.75),
                            "I"  => ("Iodine"	    ,126.9045),
                            "Te" =>	("Tellurium"	,127.60),
                            "Xe" =>	("Xenon"	    ,131.30),
                            "Cs" =>	("Cesium"	    ,132.9054),
                            "Ba" =>	("Barium"	    ,137.33),
                            "La" =>	("Lanthanum"	,138.9055),
                            "Ce" =>	("Cerium"	    ,140.12),
                            "Pr" =>	("Praseodymium"	,140.9077),
                            "Nd" =>	("Neodymium"	,144.24),
                            "Pm" =>	("Promethium"	,145),
                            "Sm" =>	("Samarium"	    ,150.4),
                            "Eu" =>	("Europium"	    ,151.96),
                            "Gd" =>	("Gadolinium"	,157.25),
                            "Tb" =>	("Terbium"	    ,158.9254),
                            "Dy" =>	("Dysprosium"	,162.50),
                            "Ho" =>	("Holmium"	    ,164.9304),
                            "Er" =>	("Erbium"	    ,167.26),
                            "Tm" =>	("Thulium"	    ,168.9342),
                            "Yb" =>	("Ytterbium"	,173.04),
                            "Lu" =>	("Lutetium"	    ,174.967),
                            "Hf" =>	("Hafnium"	    ,178.49),
                            "Ta" =>	("Tantalum"	    ,180.9479),
                            "W"  => ("Tungsten"	    ,183.85),
                            "Re" =>	("Rhenium"	    ,186.207),
                            "Os" =>	("Osmium"	    ,190.2),
                            "Ir" =>	("Iridium"	    ,192.22),
                            "Pt" =>	("Platinum"	    ,195.09),
                            "Au" =>	("Gold"	        ,196.9665),
                            "Hg" =>	("Mercury"	    ,200.59),
                            "Tl" =>	("Thallium"	    ,204.37),
                            "Pb" =>	("Lead"	        ,207.2),
                            "Bi" =>	("Bismuth"	    ,208.9804),
                            "Po" =>	("Polonium"	    ,209),
                            "At" =>	("Astatine"	    ,210),
                            "Rn" =>	("Radon"	    ,222),
                            "Fr" =>	("Francium"	    ,223),
                            "Ra" =>	("Radium"	    ,226.0254),
                            "Ac" =>	("Actinium"	    ,227.0278),
                            "Pa" =>	("Protactinium"	,231.0359),
                            "Th" =>	("Thorium"	    ,232.0381),
                            "Np" =>	("Neptunium"	,237.0482),
                            "U"  => ("Uranium"	    ,238.029),
                            "Pu" =>	("Plutonium"	,242),
                            "Am" =>	("Americium"	,243),
                            "Bk" =>	("Berkelium"	,247),
                            "Cm" =>	("Curium"	    ,247),
                            "No" =>	("Nobelium"	    ,250),
                            "Cf" =>	("Californium"	,251),
                            "Es" =>	("Einsteinium"	,252),
                            "Hs" =>	("Hassium"	    ,255),
                            "Mt" =>	("Meitnerium"	,256),
                            "Fm" =>	("Fermium"	    ,257),
                            "Md" =>	("Mendelevium"	,258),
                            "Lr" =>	("Lawrencium"	,260),
                            "Rf" =>	("Rutherfordium",261),
                            "Bh" =>	("Bohrium"	    ,262),
                            "Db" =>	("Dubnium"	    ,262),
                            "Sg" =>	("Seaborgium"	,263))
#####





## Convenience functions for mass 

atom_id(name)   = findfirst(atomic_number.==name)
id_atom(number) = atomic_number[number]

u_in_g(x)       = x * 1.66054e-24
mass_u(name::Symbol) = last(atomic_mass[String(name)])
mass_g(name::Symbol) = name |> mass_u |> u_in_g
mass_u(number::Int)  = mass_u(number |> id_atom) 
mass_g(number::Int)  = mass_g(number |> id_atom)

full_name(name) = first(atomic_mass[name])


## Plot default setup that kind-of works
#=
import PyCall.rc
rc("mathtext", fontset="cm", default="regular")
rc("font", family="monospace", size=12)
rc("text", usetex=false)
rc("xtick", direction="in", top=true)
rc("ytick", direction="in", right=true)
rc("xtick.minor", visible=true)
rc("ytick.minor", visible=true)
=#