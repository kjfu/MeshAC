using LinearAlgebra, ScatteredInterpolation, JuLIP

########## dislocation helper function ###########
function modifyA_disl(at::Atoms{Float64}, p::Real)
    xm = mat(positions(at))
    idx = findall(abs.(xm[1,:]) .>= p)
    deleteat!(at, idx)
    Xat = positions(at)
    return at, Xat
end

function idx_ref_disl(atref::Atoms{Float64}, xa_min::Float64)
    xaref = mat(positions(atref))
    I = Int64.(zero(xaref))[:]
    i = findall(xa_min-0.1 .<= xaref[:] .< -xa_min-0.1)
    I[i] .= 1
    I = reshape(I, 3, length(atref))
    idx = findall(vec(sum(I, dims=1)) .== 3)
    return idx
end
####################

function construct_ref(Lref::Real; sp=:W, defects=:SingVac, n_crack=3)
    at_ref, xcidx_ref = construct_at(Lref)
    at_max, xcidx_max = construct_at(Lref+6)
    # defect
    @show at_ref[xcidx_ref]
    @show at_max[xcidx_max]
    if defects == :SingVac
        deleteat!(at_max, xcidx_max)
        deleteat!(at_ref, xcidx_ref)
    elseif defects == :MicroCrack
        idx_max = xcidx_max-n_crack:xcidx_max+n_crack
        idx_ref = xcidx_ref-n_crack:xcidx_ref+n_crack
        # @show at_ref[idx_ref]
        # @show at_max[idx_max]
        deleteat!(at_max, idx_max)
        deleteat!(at_ref, idx_ref)
    end
    i_free = idx_ref_(at_max, at_ref)
    i_mask_zero = setdiff(1:length(at_max), i_free)
    # boundary condition
    set_pbc!(at_max, false)
    mask = fill(false, 3, length(at_max))
    mask[:, i_free] .= true
    set_mask!(at_max, mask)
    return at_max, at_ref, i_free
end

# hieracical construct A region from reference configuration
function idx_ref_(atref::Atoms{Float64}, at::Atoms{Float64})
    xaref = mat(positions(atref))
    I = Int64.(zero(xaref))[:]
    xa_min = positions(at)[1][1]
    i = findall(xa_min-0.1 .<= xaref[:] .< -xa_min-0.1)
    I[i] .= 1
    I = reshape(I, 3, length(atref))
    idx = findall(vec(sum(I, dims=1)) .== 3)
    @assert length(idx) == length(at)
    at_test = Atoms(:W, xaref[:, idx])
    @assert isapprox(mat(positions(at_test)), mat(positions(at)))
    return idx
end

# for perfecet lattice
function construct_at(L::Real)
    at = bulk(:W, cubic=true) * L
    # locate the 'center' of atomistic lattice
    Xat = positions(at)
    xcell = diag(cell(at))/2
    xcidx = findall(x->isapprox(x, xcell), Xat)
    if isnothing(xcidx)
        r = [ norm(x - xcell) for x in Xat ]
        _, xcidx = findmin(r)
        Xcidx = xcidx.I[2]
        @assert !isnothing(xcidx)
    end
    xc = Xat[xcidx]
    Xat .-= xc
    set_positions!(at, Xat)
    set_pbc!(at, false)
    return at, xcidx[1]
end

# for 3 dimentional interpolation. MATLAB implementation, very fast
# Just use for test
function fast_interp(Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2})
    write_matfile(string("./test.mat"); Xlast=Xlast, Xnext=Xnext, Ulast=Ulast)
    run(`matlab -nosplash -nodesktop -r run_matlab_test`)
    mf = MatFile(string("./test.mat")) 
    Unext = get_variable(mf, "Unext")
    close(mf)
    return Unext
end

# for 3 dimentional interpolation
# MATLAB implementation, very fast
# mxcall bug
function fast_interp(path::String, Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2})
    write_matfile(string("./data/", path, "/before.mat"); Xlast=Xlast, Xnext=Xnext, Ulast=Ulast)
    if path == "SingVac"
        run(`matlab -nosplash -nodesktop -r run_matlab_singvac`)
    elseif path == "MicroCrack"
        run(`matlab -nosplash -nodesktop -r run_matlab_microcrack`)
    else
        error("Haven't implemeted yet!")
    end
    mf = MatFile(string("./data/", path, "/after.mat")) 
    Unext = get_variable(mf, "Unext")
    close(mf)
    return Unext
end

# for 3 dimentional interpolation
# Julia implementation, very slow
function interp(Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2})
    Unext = zero(Xnext)
    for i = 1:3
        F = interpolate(Multiquadratic(), Xlast, Ulast[i,:])
        Unext[i, :] = evaluate(F, Xnext)
    end
    return Unext
end

# the final interpolation u_N -> u_int
# to test the atomistic convergence
# for 3 dimentional interpolation
function fast_interp_final(path::String, Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2}; L=15, c=0.5)
    Unext = zero(Xnext)
    Uat = load(string("./data/", path, "/L", L, ".jld"), "Uat")
    r = sqrt.(sum(abs2, Xnext, dims=1))'
    ifree = findall(r .<= c*rnn(:W)*L)
    Unext[:, ifree] .= Uat[:, ifree]
    return Unext
end

# faster interpolation: interpolation only in continuum region.
function fast_interp_ref(Xlast::Array{Float64,2}, atc, at_ref::Atoms{Float64}, Ulast::Array{Float64,2})
    Unext = zero(mat(positions(at_ref)))
    idx = atc.data["idx"]
    irefC = setdiff(1:length(at_ref), idx)
    iC = setdiff(1:size(atc.X, 2), 1:length(idx))
    Unext[:, idx] .= Ulast[:, 1:length(idx)]
    Unext[:, irefC] .= fast_interp(Xlast[:, iC], mat(positions(at_ref))[:, irefC], Ulast[:, iC])
    return Unext
end


function meshgrid(xx::Array{Float64}, yy::Array{Float64}, zz::Array{Float64})
	m, n, o = length(xx), length(yy), length(zz)
	xx = reshape(xx, 1, n, 1)
    yy = reshape(yy, m, 1, 1)
    zz = reshape(zz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
	return adjoint([vec(xx[om, :, oo]) vec(yy[:, on, oo]) vec(zz[om, on, :])])[:,:]
end

function meshgrid(xx::Array{Int}, yy::Array{Int}, zz::Array{Int})
	m, n, o = length(xx), length(yy), length(zz)
	xx = reshape(xx, 1, n, 1)
    yy = reshape(yy, m, 1, 1)
    zz = reshape(zz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
	return adjoint([vec(xx[om, :, oo]) vec(yy[:, on, oo]) vec(zz[om, on, :])])[:,:]
end

function GausWeights(ww::Array{Float64}, ii::Array{Int})
    nW = size(AA,2)
    W = zeros(nW, 1)
    for i = 1:nW
        W[i] = prod(ww[ii[:,i]])
    end
    return W
end 

function Tet_midpoint(X::Array{Float64,2})
    mpt = Array{Float64,2}(undef, 3, 6)
    for (i, j, k) ∈ zip([1,1,1,2,2,3],[2,3,4,3,4,4],[1,2,3,4,5,6])
        mpt[:,k] = 0.5*(X[:,i]+X[:,j])
    end
    mpIdx = [[1,2,3], [1,4,5], [2,4,6], [3,5,6]]
    return mpt, mpIdx
end

function Tet_circumcenter_face(X::Array{Float64,2})
    tfcc = Array{Float64,2}(undef, 3, 4)
    for (i,j,k,l) ∈ zip([1,1,1,2], [2,2,3,3], [3,4,4,4], [1,2,3,4])
        tfcc[:,l] = circumcenter_tri(X[:,[i,j,k]])
    end
    tfccIdx = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]
    return tfcc, tfccIdx
end

function circumcenter_tri(X::Array{Float64,2})
    x1, y1, z1 = X[:,1]
    x2, y2, z2 = X[:,2]
    x3, y3, z3 = X[:,3]
    a1 = y1*z2 - y1*z3 - z1*y2 + z1*y3 + y2*z3 - y3*z2
    b1 = -x1*z2 + x1*z3 + z1*x2 - z1*x3 -x2*z3 + x3*z2
    c1 = x1*y2 - x1*y3 - y1*x2 + y1*x3 + x2*y3 -x3*y2
    d1 = -x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x3*y1*z2 - x2*y3*z1 + x3*y2*z1
    a2 = 2*(x2 - x1)
    b2 = 2*(y2 - y1)
    c2 = 2*(z2 - z1)
    d2 = x1^2 + y1^2 + z1^2 - x2^2 -y2^2 - z2^2
    a3 = 2*(x3 - x1)
    b3 = 2*(y3 - y1)
    c3 = 2*(z3 - z1)
    d3 = x1^2 + y1^2 + z1^2 - x3^2 -y3^2 - z3^2

    A = [a1 b1 c1; a2 b2 c2; a3 b3 c3]
    b = [d1,d2,d3]
    return -A\b
end

function Tet_circumcenter(X::Array{Float64,2})
    r1 = sum(abs2.(X[:,1]))/2
    r2 = sum(abs2.(X[:,2]))/2
    r3 = sum(abs2.(X[:,3]))/2
    r4 = sum(abs2.(X[:,4]))/2
    A = X[:,[2,3,4]] - repeat(X[:,1], 1, 3)
    b = [r2-r1,r3-r1, r4-r1]
    return A'\b
end

function classify_tet(Tet::Array{Int64,2}, β::Vector{Float64})
    nT = size(Tet,2)
    Ttype = Vector{Int}(undef, nT)
    for i = 1:nT
        t = Tet[1:4, i]
        if all(x->x==1.0, β[t])
            Ttype[i] = 1
        elseif all(x->x==0.0, β[t])
            Ttype[i] = 0
        else
            Ttype[i] = -1
        end
    end
    return Ttype
end