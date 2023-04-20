"""
`modified from module JuLIP.Solve`

Contains a few geometry optimisation routines, for now see
the help for these:

* `minimise!`
"""
# module Solve

import Optim
import LineSearches

using Optim: OnceDifferentiable, optimize, ConjugateGradient, LBFGS
using LineSearches: BackTracking

using LinearAlgebra, SparseArrays

include("AtC.jl")
# include("AtCpro.jl")

# export minimise!

import JuLIP: minimise!

"""
`minimise!(at::AbstractAtoms)`: geometry optimisation

`at` must have a calculator and a constraint attached.

## Keyword arguments:
* `precond = :auto` : preconditioner; more below
* `gtol = 1e-6` : gradient tolerance (max-norm)
* `ftol = 1e-32` : objective tolerance
* `Optimiser = :auto`, `:auto` should always work, at least on the master
   branch of `Optim`; `:lbfgs` needs the `extraplbfgs2` branch, which is not
   yet merged. Other options might be introduced in the future.
* `verbose = 1`: 0 : no output, 1 : final, 2 : iteration and final
* `robust_energy_difference = false` : if true use Kahan summation of site energies
* `store_trace = false` : store history of energy and norm of forces
* `extended_trace = false`: also store full history of postions and forces
* `maxstep = Inf`: maximum step size, useful if initial gradient is very large
* `callback`: callback function to pass to `optimize()`, e.g. to use alternate convergence criteria

## Preconditioner

`precond` may be a valid preconditioner, e.g., `I` or `Exp(at)`, or one of
the following symbols

* `:auto` : the code will make the best choice it can with the avilable
   information
* `:exp` : will use `Exp(at)`
* `:id` : will use `I`
"""
function minimise!(config::AtC{Float64}, calc::AbstractCalculator;
                  method = :auto,
                  gtol=1e-5, ftol=1e-32,
                  verbose = 1,
                  store_trace = false,
                  extended_trace = false,
                  maxstep = Inf,
                  callback = nothing,
                  g_calls_limit = 1_000)

   # create an objective function
   obj_f = x->energy(config, x, calc)
   obj_g! = (g, x) -> copyto!(g, gradient(config, x, calc))

   precond = I
   optimiser = ConjugateGradient(linesearch = BackTracking(order=2, maxstep=maxstep))

   results = optimize( obj_f, obj_g!, get_x(config), optimiser,
                        Optim.Options( f_tol = ftol, g_tol = gtol,
                                       g_calls_limit = g_calls_limit,
                                       store_trace = store_trace,
                                       extended_trace = extended_trace,
                                       callback = callback,
                                       show_trace = (verbose > 1)) )
   # set_dofs!(at, Optim.minimizer(results))
   update!(config, Optim.minimizer(results))
   # analyse the results
   if verbose > 0
      println(results)
   end
   return results
end

function laplace_matrix(config::AtC{Float64})
   T = config.T
   Edge = [T[1,:] T[2,:]; 
         T[1,:] T[3,:];
         T[1,:] T[4,:];
         T[2,:] T[1,:];
         T[2,:] T[3,:];
         T[2,:] T[4,:];
         T[3,:] T[1,:];
         T[3,:] T[2,:];
         T[3,:] T[4,:];
         T[4,:] T[1,:];
         T[4,:] T[2,:];
         T[4,:] T[3,:];]

   Edge = unique(Edge, dims=1)

   I, J, Z = Int[], Int[], Float64[]
   for i in range(1, stop=size(Edge,1))
      ia, ja = Edge[i,:]
      append!(I, [3*(ia-1)+1, 3*(ia-1)+2, 3*(ia-1)+3])
      append!(J, [3*(ja-1)+1, 3*(ja-1)+2, 3*(ja-1)+3])
      append!(Z, -1*ones(Float64, 3, 1))
   end
   nX = 3*size(config.X,2)
   AM = sparse(I, J, Z, nX, nX)
   @assert issymmetric(AM)
   for i = range(1, stop=nX)
      AM[i,i] = -sum(AM[i,:])
   end
   return AM
end

function nsoli_condition_gradient(config::AtC{Float64}, U::Array{Float64,1}, calc::AbstractCalculator, P::SparseMatrixCSC)
   if "xfree" in keys(config.data)
		xfree = convert(Array{Int64,1}, config.data["xfree"])
	else 
		xfree = get_free!(config)
	end
	Ux = zero(config.U)
	Ux[xfree] = U
	Frc = forces(update!(config, Ux, Val{:U}()), calc, Val{:BGFC}())
   Frc = rmul!(Frc[xfree], -1.0)
   if size(P,1) > length(xfree)
      P = P[xfree, xfree]
   end
   @assert size(P,1) == length(Frc)
   return P \ Frc
end

# end
