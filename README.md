# `TSO.jl`

`TSO.jl` is a package for handling opacity + EoS tables in various formats. It is capable of reading, converting and binning monochromatic opacities on Density-Energy/Temperature grids and converting them to ouput formats, like e.g. the dispatch square gas tables. The working principle is that the code 

1. Reads the table in a supported format
2. converts to an easy, cross-platform, internal representation
3. Does all kinds of things with it
4. converts it to an output format
  
 It exports several top-level functions, that have specialized calls depending on the given table. Below is a start-to-end example on how to use the code.

_______

# Shortcuts
You can jump through the documentation by using the follwing topic shortcuts. Note that less relevant sections of the documentation may not appear here.

1. [Tabulated Spectral Opacities](#tabulated-spectral-opacities)
    1. [MARCS](#marcs)
    1. [MULTI3D](#multi3d)
    2. [AESOPUS](#aesopus)
2. [Opacity Binning](#opacity-binning)
    1. [Binning Principle](#binning-principle)
    2. [Binning Techniques](#binning-techniques)

-----------
# Tabulated Spectral Opacities

For creating 3D model atmospheres with DISPATCH the EoS and opcity table are of extreme significance. 3D RHD simulations have to -- due to the complexity of the problem -- rely on binned radiative transfer of pre-computed opacity tables. The tables that are eventually binned can have various sources, in `TSO.jl` a couple of these sources are implemented and can even be created throuch wrappers. In the following currently available formats are summarized. More information can be found in the `examples` folder, or directly in the doc strings of the corresponding functions.

## MARCS
The standard EoS + Opacity table used in DISPATCH is the `MARCS` table. A table in this format lists opacities for $\rm \sim 10^5$ wavelength points from UV to infra-red and furthermore contains information about the state of the gas.

### Creating EoS
*Relevant examples: `examples/creating_tables/create_eos.jl`*

A opacity table in the MARCS format does not contain an equation of state (EoS) with respect to the internal energy. For this purpose, one at the moment has to rely on a different code, which in our case is either Turbospectrum or MULTI3D. The following example will use Turbospectrum (TS) for this, which you can download [here](https://github.com/bertrandplez/Turbospectrum_NLTE). Please also download the Turbospectrum wrapper, which is linked on the github page as well. To aviod changing paths, put them in the same folder as `TSO.jl`. Then, git checkout the `opacity_tables` branch in TS and the wrapper. Compile TS afterwards.

If everything is setup correctly, you can run TS from within `TSO.jl`!. Create an empty project and simply specify where to find it.

```julia
using TSO

TSO.load_TS()         # Set the TS path
TSO.load_wrapper()    # Set the Wrapper path
TSO.import_wrapper()  # import the python modules of the wrapper
```

where you can specify custom paths here, if you installed it at a different location. You can then create the table, or the input files for the table, as columns using 

```julia
TSO.write_as_stagger(TColumn, RhoColumn)
```

You can follow the example `examples/creating_table/create_eos.jl` for the details. It uses the `setup = TSO.computeOpac.setup(args...)` function from the TS wrapper, to the run TS with a given chemical composition `abundances`

```julia
TSO.babsma!(setup, abundances)
```
A `bsyn!` function is also available, so it is also possible to create the entire opacity table in TS. Note however, that at the current time TS has the limitation of uniform $\rm \lambda$ steps.

### Converting Opacities
*Relevant examples: `examples/converting_tables/marcs+TSO.jl`*

After an EoS has been obtained through whatever code, and is prensent in the `TSO.SqEoS` format, it can be used to be synced with the actual opacity table from MARCS. For this purpose, there are a couple of handy functions that may be used.

```julia
using TSO
using Glob

# The EoS you created earlier
aos = @axed reload(SqEoS, abspath("../creating_table/eos_asplund07_m1_a04_v5.0.hdf5"))

# MARCS tables are split in multiple components
paths = glob("OS_table*", "opacity_tables/MARCS/asplund/Z-1.0a0.4/")

# Read raw opacity tables
mos = MARCSOpacity(paths...)

# Interpolate them to a uniform T-ρ grid in the log
m_int = uniform(mos..., new_T_size=104, new_ρ_size=104)

# Finally complement this table with the EoS
neweos, newopa, newopa_c, newopa_l, newopa_s = complement(m_int, aos, unify=false)
```

Where the `@axed` macro simply added information about the EoS axes (E or T-$\rm \rho$). `unify=false` will avoid any interpolation at this stage. As a default, only `lnEi` is taken from the provided EoS, the rest is taken from the opacity table. This means that only `lnEi` needs to be inteprolated to grid of the opacity table. You can specify this as kwargs, e.g. `complement(mos, eos; lnEi=:eos, lnRoss=:opacity)`. To remove NaN and eliminate pure outliers, you may use

```julia
set_limits!(@axed(neweos), newopa)
set_limits!(@axed(neweos), newopa_c)
set_limits!(@axed(neweos), newopa_l)
```

You can save these tables using the `save` function. It is recommended that you keep the naming convention though, so that it fits into the higher level API of some functions. It is not required though.

```julia
save(neweos,   joinpath(dname, "combined_eos.hdf5"))
save(newopa,   joinpath(dname, "combined_opacities.hdf5"))
save(newopa_c, joinpath(dname, "combined_Copacities.hdf5"))
save(newopa_l, joinpath(dname, "combined_Lopacities.hdf5"))
save(newopa_s, joinpath(dname, "combined_Sopacities.hdf5"))
```

where `C, L, S` stand for continuum, line and scattering opacties, respectively. Tables can always be reloaded using the `reload` function

```julia
eos = reload(SqEoS, "combined_eos.hdf5")
opa = reload(SqOpacity, "combined_opacities.hdf5", mmap=true)
```
where you can, in the case of opacities, specify `mmap=true` to load the memory intensive parts -- the opacities -- as memory maps to save RAM.


## MULTI3D

*Relevant examples: `MUST.jl/examples/running_multi/opacity_tables.jl`*

As a new feature, it is now also possible to create the EoS + Opacity table very conveniently in MULTI3D. For this, [MUST.jl](https://github.com/pe1995/MUST.jl) is needed, because it contains the functionality to directly run MULTI3D from within julia. Please follow the installation instructions given there. Once `MUST.jl` is installed, you need only a working MULTI3D@Dispatch installation.

From here, the idea is identical to [Turbospectrum](#creating-eos), as you create input models and then run M3D with them.

```julia
using MUST
using TSO

# load the M3D installation
MUST.@import_m3dis "path/to/Multi3D"

# specify linelists to include
linelists = String[
  "./input_multi3d/vald_2490-25540.list"
]

# Create the input for M3D, create folder if non-existing
modelatmosfolder = "input_multi3d/test_opac_table/"
models = TSO.opacityTableInput(
  MUST.@in_m3dis(modelatmosfolder),
  lnT=range(log(1.1e3), log(5.5e5); length=59) |> collect, 
  lnρ=range(log(1e-15), log(1e-3); length=59) |> collect
)
```

The opacity table can then be computed in the generic running scheme presented in the `MUST.jl` documentation. It is in principle as simple as this

```julia
opacityTable(models; 
  folder, linelist, λs, λe, δλ, 
  in_log=true, slurm=false, kwargs...) = begin
  
  MUST.whole_spectrum(
    models, 
    namelist_kwargs=(
      :model_folder=>folder,
      :linelist=>nothing,
      :absmet=>nothing,
      :linelist_params=>(:line_lists=>linelist,),
      :atom_params=>(:atom_file=>"input_multi3d/atoms/atom.h20", ),
      :spectrum_params=>(:daa=>δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>in_log),
      :atmos_params=>(
        :dims=>1, 
        :atmos_format=>"Text",
        :use_density=>true, 
        :use_ne=>false
      ),
      :m3d_params=>(
      	:n_nu=>1, 
      	:ilambd=>0,
      	:quad_scheme=>"disk_center",
      	:long_scheme=>"disk_center",
      	:make_opac_table=>true
      ),
      kwargs...
    ),
    slurm=slurm
  )
end
```

Followed by a call of this function for the chosen setup. This function as also available in `MUST.jl`.

```julia
MUST.opacityTable(
  models; 
  folder=modelatmosfolder, 
  linelist=linelists,
  λs=log(1000), λe=log(100000), δλ=0.0003, in_log=true,
  slurm=true
)
```

The output can be read, and combined together after the runs have concluded

```julia
# read the output using m3dis python routines
m3dis_models = [MUST.M3DISRun("data/$(m)") for m in models]

# collect the opacities (+EoS) from all output files
eos, opa = TSO.collect_opacity(m3dis_models)
```

They can then simply be saved using the `save` function, as mentioned above.


## AESOPUS

There is not-yet a generic, top-level API for reading various tables (although easy to add, TBD). The following example works for tables provided by the AESOPUS 2.0 opacity code. The following functions depend on the input format of the tables. The first step is to read in all different sub-tables that make up the final table

```julia
pathE = "ener_z0.0165_x0.7351_MBS22_n1458.tab"
pathK = "kross_z0.0165_x0.7351_MBS22_n1458.tab"
pathN = "nelec_z0.0165_x0.7351_MBS22_n1458.tab"
pathP = "pgas_z0.0165_x0.7351_MBS22_n1458.tab"
pathS = "sigma_z0.0165_x0.7351_MBS22_n5365.tab"
```
```julia
eosAESO = AesopusEoS(energy=pathE, kross=pathK, pg=pathP, ne=pathN)
opaAESO = AesopusOpacity(eosAESO, sigma=pathS)
```

In `eosAESO` and `opaAESO` the tables in the intermediate format are saved. The intermediate format requires that the tables are always a function of density and either internal energy or temperature. The tables are __not__ required to be uniform yet. Depending on the dimensions and approximate spacing of the arrays the code will figure that out itself and call the corresponding interpolation functions. If the data is __gridded__, the corresponding axis needs to be a 1D array, if it is furthermore evenly spaced it is called "uniform". To make a table uniform, the `uniform` function can be used, to make it gridded the `gridded` function may be used. It is futhermore possible to replace an existing EoS with a new one by using the `complement` functions. These functions can furthermore switch the energy grid from T to E and vice-versa, which can also be accomplished by the `switch_energy` function. A call like the following for example can 
  (a) switch the energy grid from T to E
  (b) make the new grid uniform, and
  (c) upsample the new energy axis to a larger number of points

```julia
eosAESO_e, opaAESO_e = uniform(eosAESO, opaAESO, upsample=5000, switch=true)
```

Keep in mind however that, depending on the size of the monochromatic opacity grid, this may create a gigantic table. It may be better to work with a T grid until after the binning. The EoS can be swapped with the command

```julia
opa_new = complement(eosAESO, eos_new, opaAESO)
```

To save the tables one may use the `save` function, where the only format currently supported is "HDF5"

```julia
save(eosAESO, "combined_eos_AESOPUS_TSO.hdf5")
```

-----------
# Opacity Binning
The binning procedure consists of 2 steps, that can be seen in the `examples/binning_opacities` folder. First, the formation height of the different wavelength need to be computed. For this, a 1D model is needed, which should be read into the generic `Model1D` type. Afterwards, the actual binning can take place, where different binning schemes can be picked.

## Binning Principle

The mathematical background of the opacity binning is described in Eitner et al. 2023. In basic terms, a large, wavelength resolved list of tabulated opacities, from whatever source, is split into a small number of wavelength bins to avoid the monochromatic solution of the radiative transfer equation during the simulation of a 3D RHD model atmosphere. This is not only a tremendous computational speed-up, but currently the only feasable way, since the RT solution is the most time consuming part in total. 

There are multiple ways opacities can be grouped together. `TSO.jl` implements a variety of different techniques, such as e.g. optical depth dependend, but wavelength independent, (MURaM) or wavelength and optical depth dependent (Stagger, Co5bold). `TSO.jl` implements those specific, but also a convenient and flexible generic binning technique that can be fine tuned easily.   

## Binning Techniques
*Relevant examples: `examples/binning_opacities/opacity_tables/bin_table.jl`*

In the given example the usage of the high-level interface is demonstrated. The example, and called functions, are designed such the tables with the right name are loaded at the right place. This is done to save RAM, because the monochromatic opacity table can be relatively expensive, even if loaded as memory map. Please note that you can of course work with the lower level functionality if this does not please you. In the following the binning is explained at a lower level. If this does not concern you, feel free to call the functions as given in the example script.

Binning opacities for DISPATCH requires three steps. First, the formation height of each point in the table in the desired average atmosphere needs to be computed. This formation height is relevant for grouping opacities that form in similar regions of the atmosphere together. This step is crutial. If the binning is not perfomed well enough, the average opacity will not be representative of the entire spectrum, and your final model atmosphere will be harmed.

```julia
using TSO

# the monochromatic table
opacities = reload(SqOpacity, "combined_opacities$(ext).hdf5", mmap=true)
aos = @axed reload(SqEoS, "combined_eos$(ext).hdf5")

# 1D model that is used for the initial condition
model = Average3D(av_path, logg=logg)

# compute the optical depth
τ_ross, τ_λ = optical_depth(aos, opacities, model)

# and use this for the formation height, i.e. -τ_ross(τ_λ=1)
d_ross, d_κ = formation_height(model, aos, opacities, τ_ross, τ_λ)

# collect into a opacity table, indicate the optical depth with a true at the end
formation_opacities = SqOpacity(d_κ, d_ross, opacities.src, opacities.λ, true)
```

You can also recompute the rosseland opacity based on the monochromatic table, and add this
to your EoS, if it comes from a different source

```julia
# compute τ_ross and overwrite it in the eos
rosseland_opacity!(aos.eos.lnRoss, aos, opacities; weights=ω_midpoint(opacities))
	
# make it consistend with the opacity table
transfer_rosseland!(aos.eos, opacities)
```

where simple midpoint integration weights for the wavelength are used to integrate the rosseland mean.
The formation height table can now be used in the binning itself. First, the style of bins can be chosen.

```julia
# load the tables, if not done already
eos = reload(SqEoS, "combined_ross_eos$(ext).hdf5")
opa = reload(SqOpacity, "combined_opacities$(ext).hdf5", mmap=true)
sopa = reload(SqOpacity, "combined_Sopacities$(ext).hdf5", mmap=true)
fopa = reload(SqOpacity, "combined_formation_opacities$(name_ext).hdf5", mmap=true)

# integration weights
weights = ω_midpoint(opacities)

# 1D model atmosphere, @optical adds τ_ross from given opacity table
model = @optical Average3D(eos, av_path, logg=logg) eos opa

# Binning style, default is generic kmeans ClusterBinning
ClusterBinning(
  TSO.KmeansBins; 
  opacities=opa, 
  formation_opacity=-log10.(fopa.κ_ross),
  Nbins=7, 
  kwargs...
)

# Sort wavelength points into specified bins    
bins = binning(bins, opacities, -log10.(fopa.κ_ross)) 
```

The last step does the sorting already. The binning of corresponding opacities and source function can then be done by simply calling the `tabulate()` function.

```julia
# scattering (sopa) can be left out
binned_opacities = tabulate(bins, weights, eos, opa, sopa, transition_model=model)
```

If no scattering opacity table is given, it is not substracted from the Planck mean and included in the opacity in the opticall thin regime (if you added it to the total opacity when creating the monochromatic table of course. If you did not, then you don't need to specify the scattering opacity here). 

If the monochromatic table is on the T-rho grid (recommended!), you can convert it to the E-rho grid (while at the same time upsamping the energy axis)

```julia
aos = @axed reload(SqEoS, "eos.hdf5")
opa = reload(SqOpacity, "binned_opacities.hdf5")

# switch from T to Eint
eosE, opaE = switch_energy(aos, opa, upsample=2048);
aosE = @axed eosE

# fill NaN values and set limits (1e±30 as default)
TSO.fill_nan!(aosE, opaE)
TSO.set_limits!(aosE, opaE)

# create a new folder
!isdir(folder_new) && mkdir(folder_new) 
 
# save the E-rho table in the Dispatch format (square_gas_mod.f90)
for_dispatch(eosE, opaE.κ, opaE.src, ones(eltype(opaE.src), size(opaE.src)...), folder_new)
```
 
-------