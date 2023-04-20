
using Pkg
# absolute path please
Pkg.activate("/Users/yswang/GitHub/JuAC/")

include("../JuliaAC/AtC.jl")
include("../JuliaAC/Solve.jl")
include("../JuliaAC/utils.jl")
include("../JuliaAC/adaptive.jl")
using SparseArrays
using Isaac
using Printf

Ra, bw, Lmsh, N = 3, 2, 15, 5

atc = AtC(Ra, bw, Lmsh, 15; meshpath="/Users/yswang/Mesher3DForSJTU/build/mesher3d");
dataPath = joinpath(pathof(JuLIP)[1:end-13], "data/")
calc = JuLIP.Potentials.FinnisSinclair(dataPath*"W-pair-Wang-2014.plt", dataPath*"W-e-dens-Wang-2014.plt");

U = zero(atc.X)
update!(atc, U, Val{:U}())

let atc = atc
    for i = 1:N

        # Solve 
        println("--------------------- Solving -----------------------------")
        P = laplace_matrix(atc)

        obj_g = x -> nsoli_condition_gradient(atc, x, calc, P)
        @time x, it_hist, ierr, x_hist = nsoli(get_x(atc), obj_g; atol=1e-3, rtol=1e-3, debug = 1, lmaxit=20, maxit=20)
        @show size(it_hist, 1)
        if ierr == 0
            @printf("---->------------------ Optimization successfull ------------------------\n")
        else
            @info("Newton solver is not executed successfully!")
        end

        # Estimate
        println("--------------------- Estimating -----------------------------")
        ∇U, _, J = gradient(atc)
        Est = [norm(∇U[:,:,i] - I,2) for i in 1:size(∇U, 3)] 

        # Mark
        println("--------------------- Marking -----------------------------")
        Xn, TIdx = Adaptive.mark(atc, Est)

        # Refine
        println("--------------------- Refining -----------------------------")
        atcnew = Adaptive.refine!(atc, Xn, TIdx; meshpath="/Users/yswang/Mesher3DForSJTU/build/mesher3d") # TODO: check if atcnew == atc

        Unew = zero(atcnew.X);
        update!(atcnew, Unew, Val{:U}())

        atc = deepcopy(atcnew)

    end
end

