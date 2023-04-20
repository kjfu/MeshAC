using Printf

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
