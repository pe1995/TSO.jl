
"""Interpolate the EoS + Opacity tables based on missing values in lnEi."""
function smooth!(eos::E, opacities_list::NTuple{N,O}; along=:T) where {N, E<:EoSTable, O<:OpacityTable}
    # smooth the eos first
    m = smooth!(eos, return_missing=true, along=along)

    # return if they are all present, nothing to smooth here
    all(.!m) && return

    # check if we have an energy axis
    eaxis = ndims(eos.lnT) > 1

    for opacities in opacities_list
        # interpolate all quantities of the opacities
        # Note that, if the energy is the axis we still have to interpolate in T!
        if (along ==:T) 
            Taxis = eaxis ? zeros(size(eos.lnT, 1)) : zeros(length(eos.lnT))
            mask  = trues(length(Taxis))
            nmask = similar(mask)
            val   = similar(Taxis)

            for i in eachindex(eos.lnRho)
                # interpolate along T
                mask  .= m[:, i]
                nmask .= .!m[:, i]
                Taxis .= eaxis ? eos.lnT[:, i] : eos.lnT

                nodes = (view(Taxis, nmask),)

                val .= log.(opacities.κ_ross[:, i])
                opacities.κ_ross[mask, i] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Taxis, mask)))

                for j in eachindex(opacities.λ)
                    val .= log.(opacities.κ[:, i, j])
                    opacities.κ[  mask, i, j] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Taxis, mask)))

                    val .= log.(opacities.src[:, i, j])
                    opacities.src[mask, i, j] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Taxis, mask)))
                end

            end

        elseif (along==:Rho) | (along==:rho) | (along==:ρ) 
            eaxis && error("You can only integrate along Rho if T is the other axis (because we have to integrate along T anyways, which then may be unstructured). Choose :T.")
            
            Raxis = eos.lnRho
            mask  = trues(length(Raxis))
            nmask = similar(mask)
            val   = similar(Raxis)


            for i in eachindex(eos.lnT)
                # interpolate along T
                mask  .= m[i, :]
                nmask .= .!m[i, :]

                nodes = (view(Raxis, nmask),)

                val .= log.(opacities.κ_ross[i, :])
                opacities.κ_ross[i, mask] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Raxis, mask)))

                for j in eachindex(opacities.λ)
                    val .= log.(opacities.κ[i, :, j])
                    opacities.κ[  i, mask, j] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Raxis, mask)))

                    val .= log.(opacities.src[i, :, j])
                    opacities.src[i, mask, j] .= exp.(extrapolate(interpolate(nodes, view(val, nmask), Gridded(Linear())), Line()).(view(Raxis, mask)))
                end
            end
        else
            error("Given axis not implemented. Choose T or Rho")
        end
    end

    nothing
end

"""Interpolate the EoS tables based on missing values in lnEi."""
function smooth!(eos::E; along=:T, return_missing=false) where {E<:EoSTable}
    # missing spots in the eos
    m = missing_spots(eos)

    # return if they are all present, nothing to smooth here
    all(.!m) && return

    # check if we have an energy axis
    eaxis = ndims(eos.lnT) > 1

    # interpolate all quantities of the eos
    # Note that, if the energy is the axis we still have to interpolate in T!
    if (along ==:T) 
        Taxis = eaxis ? zeros(size(eos.lnT, 1)) : zeros(length(eos.lnT))
        mask  = trues(length(Taxis))
        nmask = similar(mask)

        for i in eachindex(eos.lnRho)
            # interpolate along T
            mask  .= m[:, i]
            nmask .= .!m[:, i]
            Taxis .= eaxis ? eos.lnT[:, i] : eos.lnT

            nodes = (view(Taxis, nmask),)
            eos.lnEi[  mask, i] .= extrapolate(interpolate(nodes, view(eos.lnEi,   nmask, i), Gridded(Linear())), Line()).(view(Taxis, mask))
            eos.lnPg[  mask, i] .= extrapolate(interpolate(nodes, view(eos.lnPg,   nmask, i), Gridded(Linear())), Line()).(view(Taxis, mask)) 
            eos.lnNe[  mask, i] .= extrapolate(interpolate(nodes, view(eos.lnNe,   nmask, i), Gridded(Linear())), Line()).(view(Taxis, mask)) 
            eos.lnRoss[mask, i] .= extrapolate(interpolate(nodes, view(eos.lnRoss, nmask, i), Gridded(Linear())), Line()).(view(Taxis, mask))  
        end
    elseif (along==:Rho) | (along==:rho) | (along==:ρ) 
        eaxis && error("You can only integrate along Rho if T is the other axis (because we have to integrate along T anyways, which then may be unstructured). Choose :T.")
        
        Raxis = eos.lnRho
        mask  = trues(length(Raxis))
        nmask = similar(mask)

        for i in eachindex(eos.lnT)
            # interpolate along T
            mask  .= m[i, :]
            nmask .= .!m[i, :]

            nodes = (view(Raxis, nmask),)
            eos.lnEi[  i, mask] .= extrapolate(interpolate(nodes, view(eos.lnEi,   i, nmask), Gridded(Linear())), Line()).(view(Raxis, mask))
            eos.lnPg[  i, mask] .= extrapolate(interpolate(nodes, view(eos.lnPg,   i, nmask), Gridded(Linear())), Line()).(view(Raxis, mask)) 
            eos.lnNe[  i, mask] .= extrapolate(interpolate(nodes, view(eos.lnNe,   i, nmask), Gridded(Linear())), Line()).(view(Raxis, mask)) 
            eos.lnRoss[i, mask] .= extrapolate(interpolate(nodes, view(eos.lnRoss, i, nmask), Gridded(Linear())), Line()).(view(Raxis, mask))  
        end
    else
        error("Given axis not implemented. Choose T or Rho")
    end

    return if return_missing
        m
    else
        nothing
    end
end

smooth!(eos::E, opacities::O, args...; kwargs...) where {E<:EoSTable, O<:OpacityTable} = smooth!(eos, (opacities,), args...; kwargs...)

"""Get the missing spots in the EoS."""
function missing_spots(eos::E) where {E<:EoSTable}
    m = falses(size(eos.lnEi)...)

    for j in axes(m, 2)
        for i in axes(m,1)
            m[i, j] = isnan(eos.lnEi[i, j])
        end
    end

    m
end

smooth(eos::E, args...; kwargs...) where {E<:EoSTable} = begin
    eos_copy = deepcopy(eos)
    smooth!(eos_copy, args...; kwargs...)
    
    eos_copy
end

smooth(eos::E, opacities::O, args...; kwargs...) where {E<:EoSTable, O<:OpacityTable} = begin
    eos_copy = deepcopy(eos)
    opa_copy = deepcopy(opacities)
    smooth!(eos_copy, opa_copy, args...; kwargs...)
    
    eos_copy,opa_copy
end