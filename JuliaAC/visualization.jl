"""
## module Visualization

## Description:
1. plot dealII geometry i.e. cubes
""" 
module Visualization

using PyPlot
using JuLIP

function cube_dealII(X::Array{Float64,2}, T::Array{Int64,2})
    Ii = [1 2 3 4 1 2 3 4 5 6 7 8; 
          2 3 4 1 5 6 7 8 6 7 8 5]
    for iT = 1:size(T, 2)
        tt = T[:, iT]
        for iI = 1:12
            i = Ii[:, iI]
            xx = X[:, tt[i]]
            plot3D(xx[1,:], xx[2,:], xx[3,:], "bo-", linewidth=0.7, markersize=4)
        end
    end
    xticks([])
    yticks([])
    zticks([])
end

function scatter_atom(at::Atoms)
    X = mat(at.X)
    scatter3D(X[1,:], X[2,:], X[3,:])
end

function tetrahedron(X::Array{Float64,2})
    Ii = [1 1 1 2 2 3; 
          2 3 4 3 4 4]

    for i = 1:size(Ii,2)
        ii = Ii[:,i]
        xx = X[:,ii]
        plot3D(xx[1,:], xx[2,:], xx[3,:], "bo-", linewidth=0.7, markersize=4)
    end

    xticks([])
    yticks([])
    zticks([])
end

end # End module