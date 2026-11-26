<img src="tso_logo.png" height=100% width=100%/> 

`TSO.jl` is a package for handling EoS and opacity tables. It reads monochromatic opacity tables, bins them for `DISPATCH`, and provides fast lookups of thermodynamic quantities and opacities on density-temperature or density-internal energy grids. It is complementary to the `MUST.jl` package, which handles the 3D models.

The tables are created and binned with `m3dis.jl` in the `stellar_atmospheres` experiment (modes `opacity-table` and `opacity-binning`); see the M3DIS README for that. This page gives an overview of how to work with the resulting tables directly.

_______

# Shortcuts

1. [Table Folders](#table-folders)
2. [Loading Tables](#loading-tables)
3. [Looking up Values](#looking-up-values)
4. [Binning](#binning)

-----------
# Table Folders

Every table is a folder with one `HDF5` file that holds the whole table, together with everything else that belongs to it:

```
z_0.0_alpha_0.0_vmic_2.0_v4.6/
├── z_0.0_alpha_0.0_vmic_2.0_v4.6.hdf5   monochromatic table (named after the folder)
└── table_info.toml                      composition and settings (also stored in the table file)

V51MB11S_M3DIS_t05780_g+4.60_.../
├── binned_eos.hdf5                      binned table
├── *.dat, tabparam*.in                  the same table in the format read by DISPATCH
└── table_info.toml, <config>.toml, ...  settings of the binning, 1D models, figures
```

The table file holds the table on one or two grids, the density-temperature grid (`temperature`) and the density-internal energy grid used by `DISPATCH` (`energy`). Each grid has its axes and up to three parts, every part holding any number of quantities:

```
table.hdf5
├── temperature/
│   ├── lnRho, lnT     axes
│   ├── eos/           lnEi, lnPg, lnNe, lnRoss, lnK500 (500nm opacity), ...
│   ├── opacity/       κ, src, λ (wavelength or bin), ...
│   └── scattering/    κ, ... (monochromatic tables)
└── energy/            the same on the internal energy grid (binned tables)
```

The Rosseland opacity is only stored in the EoS (`lnRoss`), next to the opacity at 500 nm (`lnK500`).

The composition, version and all settings of the table are stored as metadata of the file. The file is recognised by its metadata, not by its name. Every quantity has its unit as attribute (`unit`), in Julia `units(tab)`:

| quantity | monochromatic | binned |
|---|---|---|
| `κ` | cm² g⁻¹ | cm⁻¹ (includes the density) |
| `src` | erg s⁻¹ cm⁻² sr⁻¹ cm⁻¹ (per cm of wavelength) | erg s⁻¹ cm⁻² sr⁻¹ (integrated over the bin) |
| `λ` | Å | bin index |
| `lnRoss`, `lnK500` | ln(cm² g⁻¹) | ln(cm² g⁻¹) |
| `lnRho`, `lnT`, `lnEi`, `lnPg`, `lnNe` | ln of g cm⁻³, K, erg g⁻¹, dyn cm⁻², cm⁻³ | same |

Tables in the old layout, with one file per quantity (`combined_eos_*.hdf5`, `eos_T.hdf5`, `binned_opacities.hdf5`, ...), can still be read. To combine them into one file use

```julia
migrateTable(folder)               # writes the new file next to the old ones
migrateTable(folder; delete=true)  # and removes the old files after checking the new one
```

-----------
# Loading Tables

Tables are loaded from their folder (or file) as `EquationOfState`, by default on the temperature grid:

```julia
using TSO

tab = EquationOfState("path/to/table")                  # temperature grid
tabE = EquationOfState("path/to/table"; grid=:energy)   # internal energy grid
availableGrids("path/to/table")                         # [:temperature, :energy]

tab = EquationOfState("path/to/table"; parts=(:eos,))   # only the EoS
tab = EquationOfState("path/to/table"; mmap=false)      # opacities in memory instead of memory-mapped (read-only)
```

All quantities and axes are accessed by name, and any quantity can be added:

```julia
tab[:lnPg]      # or tab.lnPg; [lnT, lnRho]
tab[:κ]         # opacity [lnT, lnRho, λ or bin]
tab.λ           # wavelength (or bin index)
tab.meta        # kind, path and composition of the table

tab[:my_quantity] = zeros(size(tab)...)
```

The composition of a table is stored with it:

```julia
info = TableInfo("path/to/table")
show_composition(info)    # "[Fe/H]=-5.0, [α/Fe]=0.4, with: [C/Fe]=3.0,[N/Fe]=3.0"
```

Tables can be saved again with `TSO.save(tab, "path.hdf5")`, or both grids at once with `TSO.save("path.hdf5", tab, tabE)`.

-----------
# Looking up Values

Use `sample` with the coordinates given as keywords. Coordinates can be given linearly (`rho`, `T`, `Ei`, `Pg`) or as natural logarithm (`lnRho`, `lnT`, `lnEi`, `lnPg`), as numbers or arrays:

```julia
# one quantity
lnNe = sample(tab, :lnNe; rho=1e-7, T=6000.0)

# several quantities at once, e.g. for a whole 3D cube
lnPg, lnRoss = sample(tab, :lnPg, :lnRoss; rho=b.data[:d], T=b.data[:T])

# opacities and source functions, where k is the wavelength or bin index
κ, S = sample(tab, :κ, :src; rho=1e-7, T=6000.0, k=1)

# the scattering opacity
σ = sample(scattering(tab), :κ; rho=1e-7, T=6000.0, k=1)
```

Coordinates that are not the axes of the table are inverted, e.g. `Ei` on the temperature grid or `Pg` instead of the density.

-----------
# Binning

The opacities of a monochromatic table are binned with a bin assignment (one bin index per wavelength, or a matrix of bin weights):

```julia
binned = binning(assignment, tab.opacity.weights, tab; logg=4.44)
binned[:κ], binned[:src]   # binned opacity and source function

for_dispatch_T(binned, binned[:κ], binned[:src], ones(size(binned[:src])), "path/to/folder")
```

In practice this, together with finding the bins, is done by `m3dis.jl opacity-binning`.

-------------------------

*Claude Opus 5.5 has been used to improve text quality.*
