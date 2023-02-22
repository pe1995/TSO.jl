using Pkg; Pkg.activate("."); 
using TSO
using PyPlot
using Glob
using Serialization

begin
    paths = glob("OS_table*", "OPAC-for-3D/Z0.0a0.0");
    mos   = TSO.MARCSOpacity(paths...);

    m_int = TSO.square(mos...);

    serialize("m_int", m_int)
end
