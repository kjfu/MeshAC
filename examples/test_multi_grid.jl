using JuLIP
using LinearAlgebra
using ProgressMeter

include("../JuliaAC/FIO.jl")

Ra, bw, Rbuf = 3, 2, 2
Lmsh = 15
@assert Ra + bw + Rbuf < Lmsh

# atomistic lattice
at = bulk(:W, cubic=true)*2*(Ra+bw+Rbuf)
set_pbc!(at, false)

r0 = rnn(:W)
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
Lat = xcell[1] 

Ra *= r0
bw *= r0
Rbuf *= r0
Lmsh *= r0
cb = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1]
Xcb = Lmsh .* cb
Xcb = Xcb .+ mat(xc)

X = hcat(mat(Xat), Xcb)

Xtype = zeros(Int64, length(at))
append!(Xtype, ones(Int64, 8))

fn = "../FIO/cp.mesh"
ACFIO.write_mesh(fn, X, Xtype)

## call `mesher3d` to build coupled mesh
ofn = "../FIO/out3d.mesh"
run(`/home/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d -s 10 $fn -o $ofn`) # ML: device dependent!

X, T = ACFIO.read_mesh(ofn)
iBdry = findall(x->x==2.0, X[4,:]) 

U = rand(3, size(X,2))*1e-2
# U = zeros(3, size(X,2))
wat = bulk(:W, cubic = true)
V0=det(cell(wat))/2
DAtC = Dict{String, Any}()
include("../JuliaAC/AtC.jl")

∇U, volT, J = gradient(T, X, U)

atc = AtC{eltype(X)}(at, V0, X[1:3, :], U, ∇U, T[1:4, :], Ra, bw, J, wat, iBdry, DAtC)
atc.data["xc"] = copy(mat(xc))
atc.data["volT"] = copy(volT)
# calculator
# eam = EAM("../potentials/w_eam4.fs"; s = 1e-8)
dataPath = joinpath(pathof(JuLIP)[1:end-13], "data/")
calc = JuLIP.Potentials.FinnisSinclair(dataPath*"W-pair-Wang-2014.plt", dataPath*"W-e-dens-Wang-2014.plt")

F(config::AtC{Float64}) = energy(config, calc, Val{:BGFC}(); bfcn=:affine)   # 
dF(config::AtC{Float64}) = forces(config, calc, Val{:BGFC}(); bfcn=:affine)

include("../JuliaAC/Testing.jl")
dEh = cfdtest(F, dF, atc; verbose=true)

