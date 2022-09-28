"""Get the rosseland optical depth where the monochromatic optical depth is unity."""
function formation_height(o::RegularOpacityTable)
    @assert o.optical_depth

    z = size(o.κ, 2) # The shape of the opacity is (lambda, eos_axes...) so the second is the height scale

    h = zeros(size(o.κ, 1))
    n0::Int = 0
    @inbounds for i in eachindex(h)
        n0 = findfirst(j->o.κ[i, j]>=1.0, 1:z) - 1
        h[i] = o.κ_ross[i, n0]
    end

    h
end