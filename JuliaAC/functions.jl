include("AtC.jl")
include("FIO.jl")

using ScatteredInterpolation

function construct_geom(Ra::Int64, bw::Int64; 
                        Rbuf=2, Lmsh=15, sp=:W, r0=rnn(:W))
    
    ## construct atomistic region
    at = bulk(sp, cubic=true)*2*(Ra+bw+Rbuf)
    set_pbc!(at, false)
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
    # Lat = xcell[1] # YS: no need?

    ## construct continuum information
    Ra *= r0
    bw *= r0
    Rbuf *= r0
    Lmsh *= r0
    cb = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1]
    Xcb = Lmsh .* cb
    X = hcat(mat(Xat), Xcb)
    Xtype = zeros(Int64, length(at))
    append!(Xtype, ones(Int64, 8))
    fn = "../FIO/cp.mesh"
    ACFIO.write_mesh(fn, X, Xtype)
    # call `mesher3d` to build coupled mesh
    ofn = "../FIO/out3d.mesh"
    run(`/home/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d -s 10 $fn -o $ofn`) # ML: device dependent!
    X, T = ACFIO.read_mesh(ofn)
    iBdry = findall(x->x==2.0, X[4,:])

    ## alternative data
    data = Dict{String, Real}()
    data["xc"] = xc
    data["Ra"] = Ra
    data["bw"] = bw
    data["iBdry"] = iBdry

    return at, X, T, data

end


dataPath = joinpath(pathof(JuLIP)[1:end-13], "data/")
calc = JuLIP.Potentials.FinnisSinclair(dataPath*"W-pair-Wang-2014.plt", dataPath*"W-e-dens-Wang-2014.plt")

function energy_force_BGFC(at::Atoms{Float64}, X::Array{Float64,2}, T::Array{Int64,2}, data::Dict{String,Real}, U::Array{Float64,2}, calc::AbstractCalculator; 
                           wat=bulk(:W, cubic=true), V0=15.777248000000004)

    DAtC = Dict{String, Any}()
    ∇U, volT, J = gradient(T, X, U)

    atc = AtC{eltype(X)}(at, V0, X[1:3, :], U, ∇U, T[1:4, :], data["Ra"], data["bw"], J, wat, data["iBdry"], DAtC)
    atc.data["xc"] = data["xc"]
    atc.data["volT"] = copy(volT)

    E = energy(atc, calc, Val{:BGFC}(); bfcn=:affine)   # 
    F = forces(atc, calc, Val{:BGFC}(); bfcn=:affine)
    return E, F
end

function scatter_interp(Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2})
    F = interpolate(Multiquadratic(), Xlast, Ulast)
    Unext = evaluate(F, Xnext)
    return Unext
end



    