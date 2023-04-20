"""
## module GeomDealII

### Summary

This module contains the geometry operator of mesh constructed by deallII
1. cube_neigs
"""
module GeomDealII

using StaticArrays

"""
## function ref_cube_neigs()
### Return connection list (Dict) of node within the reference cube of deallII. Compute once and for all
"""
function ref_cube_neigs()
    I = [1 2 3 4 1 2 3 4 5 6 7 8 2 3 4 1 5 6 7 8 6 7 8 5]

    refNeigs = Dict{Int64, SVector{3, Int64}}()

    for i = 1:8
        nei = findall(x -> x==i, I)
        nei = map(x->isodd(x) ? x+1 : x-1, nei)
        nei = I[nei]
        refNeigs[i] = SVector{3, Int64}(nei)
    end
    return refNeigs
end

"""
## function cube_nodal_stress
### return deformation matrix of selected node of selected cube computed with connecting nodes.
"""
function cube_nodal_stress(refNeigs::Dict{Int64, SVector{3, Int64}}, 
                            X::Array{Float64, 2}, T::Array{Int64, 2}, U::Array{Float64, 2},
                            Tidx::Int64, Nidx::Int64)
    t = T[:,Tidx]
    xx = X[:, [t[Nidx]; collect(t[refNeigs[Nidx]])]]
    uu = U[:, [t[Nidx]; collect(t[refNeigs[Nidx]])]]
    J = zeros(Float64, 3, 3)
    Du = zeros(Float64, 3, 3)

    for j = 1:3
        J[:,j] = xx[:, j+1] - xx[:, 1]
        Du[:,j] = uu[:, j+1] - uu[:, 1]
    end
    return Du/J
end

function cube_stress(X::Array{Float64, 2}, T::Array{Int64, 2}, U::Array{Float64, 2}, 
                    Tidx::Int64, refNeigs::Dict{Int64, SVector{3, Int64}})
    Du = zeros(Float64, 3, 3)
    for i = 1:8
        Du += cube_nodal_stress(refNeigs, X, T, U, Tidx, i)
    end
    return Du/8
end

end #module