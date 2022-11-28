using Pkg; Pkg.activate(".")
using TSO

function combine(path1, extension1, path2, extension2, parent_folder=dirname(path1))

    e1  = TSO.reload(TSO.RegularEoSTable, joinpath(path1, "unified_eos_$(extension1)_step2.hdf5"), mmap=false)
    e2  = TSO.reload(TSO.RegularEoSTable, joinpath(path2, "unified_eos_$(extension2)_step2.hdf5"), mmap=false)
    
    # General opacities
    o1  = TSO.reload(TSO.RegularOpacityTable, joinpath(path1, "unified_opacities_$(extension1)_step2.hdf5"),  mmap=false)
    o2  = TSO.reload(TSO.RegularOpacityTable, joinpath(path2, "unified_opacities_$(extension2)_step2.hdf5"),  mmap=false)
    e, o = TSO.combine_opacities(e1, o1, e2, o2)
    TSO.save(e, joinpath(parent_folder, "combined_eos.hdf5"))
    TSO.save(o, joinpath(parent_folder, "combined_opacities.hdf5"))

    @info "Opacities 1/4"

    # Cont opacities
    o1  = TSO.reload(TSO.RegularOpacityTable, joinpath(path1, "unified_Copacities_$(extension1)_step2.hdf5"),  mmap=false)
    o2  = TSO.reload(TSO.RegularOpacityTable, joinpath(path2, "unified_Copacities_$(extension2)_step2.hdf5"),  mmap=false)
    e, o = TSO.combine_opacities(e1, o1, e2, o2)
    TSO.save(o, joinpath(parent_folder, "combined_Copacities.hdf5"))

    @info "Opacities 2/4"

    # Line opacities
    o1  = TSO.reload(TSO.RegularOpacityTable, joinpath(path1, "unified_Lopacities_$(extension1)_step2.hdf5"),  mmap=false)
    o2  = TSO.reload(TSO.RegularOpacityTable, joinpath(path2, "unified_Lopacities_$(extension2)_step2.hdf5"),  mmap=false)
    e, o = TSO.combine_opacities(e1, o1, e2, o2)
    TSO.save(o, joinpath(parent_folder, "combined_Lopacities.hdf5"))

    @info "Opacities 3/4"

    # Scat opacities
    o1  = TSO.reload(TSO.RegularOpacityTable, joinpath(path1, "unified_Sopacities_$(extension1)_step2.hdf5"),  mmap=false)
    o2  = TSO.reload(TSO.RegularOpacityTable, joinpath(path2, "unified_Sopacities_$(extension2)_step2.hdf5"),  mmap=false)
    e, o = TSO.combine_opacities(e1, o1, e2, o2)
    TSO.save(o, joinpath(parent_folder, "combined_Sopacities.hdf5"))
    
    @info "Opacities 4/4"
    ;
end

IR = "IR_Magg_v1.1_opacities"
UV = "UV_Magg_v1.1_opacities"
parent_folder="tables"

combine(IR, "IR_Magg_v1.1", UV, "UV_Magg_v1.1", parent_folder)