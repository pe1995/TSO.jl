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

# constants (In agreement with Tabgen)
const EiExtra    = 5.0            # extra eV per H atom to add
const KBoltzmann = 1.380658E-16   # Boltzman's cst. [erg/K]
const CLight     = 2.99792458E+10 # Speed of light [cm/s]
const HPlanck    = 6.6260755E-27  # Planck's constant [erg s]
const ev_to_erg  = 1.60218e-12    # conversion
const HIonPot    = 13.595         # hydrogen ioniz. potential [eV]
const twohc2     = 2.0e0 *HPlanck*CLight^2
const hc_k       = HPlanck*CLight/KBoltzmann
const aa_to_cm   = 1.0e-8