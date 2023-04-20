

"""
## module FIO

## Description:
File I/O stream controler
1. read geometrical contents X, T from .msh file exported by dealII
""" 


module ACFIO
using DelimitedFiles
using Printf
# read mesh information in 3D
function read_dealii_3D(fn::AbstractString)
    open(fn, "r") do io
        read_msh_dealii_3D(io)
    end
end
    
function read_msh_dealii_3D(io)
    dim = 3
    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)
    X = zeros(dim+1, NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        X[:,i] = d
    end

    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)
    T = Array{Int64, 2}(zeros(dim+10, NV))
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Int64)
        T[:,i] = d
    end
    return X[2:4, :], T[6:end, :]
end

"""
    read_mesh(filename) -> mesh::Mesh
"""
function read_mesh(fn::AbstractString)
    open(fn, "r") do io
        read_mesh(io)
    end
end

"""
    read_mesh(iostream) -> mesh::Mesh

Reads the mesh nodes, edges and elements stored in the input .mesh file

Returns an object `mesh` of type `Mesh`, comprising both vector arrays.
"""
function read_mesh(io)

    thisLine = io |> readline |> strip
    while thisLine != "Dimension"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    dim = parse(Int, thisLine)
#   println("Dimension = ", dim)
    if dim != 3
        error("Only 3 dimensional problem are considered so far!")
    end

	thisLine = io |> readline |> strip
    while thisLine != "Vertices"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)

#    P = SVector{dim+1,Float64}
    Vex = Array{Float64}(undef,dim+1,NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        Vex[:,i] = d
    end

    thisLine = io |> readline |> strip
    while thisLine != "Tetrahedra"
        thisLine = io |> readline |> strip
    end

	thisLine = io |> readline |> strip
    NT = parse(Int, thisLine)

#    P = SVector{dim+2,Int64}
    Tet = Array{Int64}(undef,dim+2,NT)
    for i in 1:NT
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Int64)
        Tet[:,i] = d
    end
    return (Vex,Tet)
end

function write_mesh(filename, X::Array{Float64,2}, Xtype::Vector{Int64})
	@assert size(X, 2) == length(Xtype)
    f = open(filename, "w")
	@printf(f, "MeshVersionFormatted 2\n\n\nDimension 3\n\n\n");
	NV = size(X, 2)
	@printf(f, "Vertices\n%d\n", NV);
	for i = 1:NV
		@printf(f, "%f %f %f %d\n", X[1,i], X[2,i], X[3,i], Xtype[i])
	end
	@printf(f, "\n\nEnd\n")
    close(f)
end

end