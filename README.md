<a href="url"><img src="tso_logo.png" height=100% width=100%/></a> 

`TSO.jl` is a package for handling EoS and opacity tables. It reads monochromatic opacity tables, bins them for `DISPATCH`, and provides fast lookups of thermodynamic quantities and opacities on density-temperature or density-internal energy grids. It is complementary to the `MUST.jl` package, which handles the 3D models.

The tables are created and binned with `m3dis.jl` in the `stellar_atmospheres` experiment (modes `opacity-table` and `opacity-binning`); see the M3DIS README for that. This page gives an overview of how to work with the resulting tables directly.

_______

# Shortcuts

1. [Table Folders](#table-folders)
2. [Loading Tables](#loading-tables)
3. [Looking up Values](#looking-up-values)

-----------
# Table Folders

There are two kinds of table folders.

A __monochromatic table__ (from `m3dis.jl opacity-table`) contains the EoS and the opacities for every wavelength point on a density-temperature grid:
  - `combined_eos_<id>.hdf5`: the EoS, including the Rosseland opacity
  - `combined_eos500_<id>.hdf5`: the same with the 500nm opacity instead
  - `combined_opacities_<id>.hdf5`: the opacity and source function for every wavelength
  - `combined_sopacities_<id>.hdf5`: the scattering opacity
  - `table_info.toml`: the chemical composition and version of the table

A __binned table__ (from `m3dis.jl opacity-binning`) contains the binned opacities, both on the density-temperature grid and on the density-internal energy grid that `DISPATCH` uses:
  - `eos_T.hdf5`, `binned_opacities_T.hdf5`: EoS and binned opacities on the temperature grid
  - `eos.hdf5`, `binned_opacities.hdf5`: the same on the internal energy grid
  - `eos_T500.hdf5`: the EoS with the 500nm opacity
  - `*.dat`, `tabparam*.in`: the same tables in the format that is read by `DISPATCH`
  - `table_info.toml`: the composition and all settings of the binning

-----------
# Loading Tables

Tables are loaded with `reload`. The opacities can be loaded as memory maps to save RAM, which is recommended for the large monochromatic tables:

```julia
using TSO

eos = TSO.reload(TSO.SqEoS, "path/to/eos_T.hdf5")
opa = TSO.reload(TSO.SqOpacity, "path/to/binned_opacities_T.hdf5", mmap=true)

opa.κ      # opacity    [ρ, T or E, λ or bin]
opa.src    # source function
opa.λ      # wavelength (or bin index)
```

The composition of a table is stored next to it:

```julia
info = TSO.TableInfo("path/to/table/folder")
TSO.show_composition(info)    # "[Fe/H]=-5.0, [α/Fe]=0.4, with: [C/Fe]=3.0,[N/Fe]=3.0"
```

Tables can be saved again with `TSO.save(eos, "path.hdf5")`.

-----------
# Looking up Values

For lookups, wrap the tables with `extended`, and use `sample` with the coordinates given as keywords. Coordinates can be given linearly (`rho`, `T`, `Ei`) or as natural logarithm (`lnRho`, `lnT`, `lnEi`), as numbers or arrays:

```julia
eos = TSO.reload(TSO.SqEoS, "path/to/eos_T.hdf5") |> TSO.extended

# one quantity
lnNe = TSO.sample(eos, :lnNe; rho=1e-7, T=6000.0)

# several quantities at once, e.g. for a whole 3D cube
lnPg, lnRoss = TSO.sample(eos, :lnPg, :lnRoss; rho=b.data[:d], T=b.data[:T])
```

Opacities and source functions are looked up together with the EoS, where `k` is the wavelength or bin index:

```julia
opa = TSO.reload(TSO.SqOpacity, "path/to/binned_opacities_T.hdf5") |> TSO.extended
κ, S = TSO.sample(eos, opa, :κ, :src; rho=1e-7, T=6000.0, k=1)
```

On a table on the internal energy grid (`eos.hdf5`), use `Ei` instead of `T`. `TSO.is_internal_energy(eos.eos)` tells you which grid a table is on.
