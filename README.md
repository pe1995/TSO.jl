# TSO.jl

TSO.jl is a package for handling opacity + EoS tables in various formats. It is capable of reading, converting and binning monochromatic opacities on Density-Energy/Temperature grids and converting them to ouput formats, like e.g. the dispatch square gas tables. The working principle is that the code 
  (a) Reads the table in a supported format
  (b) converts to an easy, cross-platform, internal representation
  (c) Does all kinds of things with it
  (d) converts it to an output format
  
 It exports several top-level functions, that have specialized calls depending on the given table. Below is a start-to-end example on how to use the code.


## Start-to-End example: AESOPUS 2.0 opacities + EoS
### Reading and Converting 

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

The table can later be reloaded using the `reload` function, e.g.

```julia
eos = reload(SqEoS, "combined_eos.hdf5")
opa = reload(SqOpacity, "combined_opacities.hdf5")
```


### Binning

The binning procedure consists of 2 steps, that can be seen in the examples folders. First, the formation height of the different wavelength need to be computed. For this, a 1D model is needed, which should be read into the generic `Model1D` type. Afterwards, the actual binning can take place, where different binning schemes can be picked.

