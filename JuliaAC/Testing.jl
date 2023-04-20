using Printf
using JuLIP
using LinearAlgebra
using ProgressMeter

include("AtC.jl")

function cfdtest(F::Function, dF::Function, config::AtC{Float64}; verbose=true)
   H = Float64[]
   errors = Float64[]
   Vex = config.X
   Tet = config.T
   # x = zeros(Float64, 3, size(Vex,2))
   x = deepcopy(config.U)
   # x = config.U
   # update!(config, x, Val{:U}())
   E = F(config)
   dE = dF(config)
   lenp = 4
   # x = zeros(Float64, 3, size(Vex,2))
   iBdry = config.iBdry
   # loop through finite-difference step-lengths
   ldE = size(dE, 2)
   prg = Progress(3*(lenp-1)*ldE)
   dEh = copy(dE)
   for p = 2:lenp
      h = 0.1^(p)
      push!(H, h)
      for n = 1:size(dE,2)
           if n ∈ iBdry
               continue
           end
           for j = 1:3
               x[j,n] += h
               # config.U = x
               update!(config, x, Val{:U}())
               Eh = F(config)
               # dEh[j,n] = (Eh - E) / h
               dEh[j,n] = (E - Eh) / h
               x[j,n] -= h
               update!(config, x, Val{:U}())
               ProgressMeter.update!(prg,3*(p-2)*ldE + 3*(n-1)+j)
           end
      end
      push!(errors, norm(dE - dEh, Inf))
   end
   @printf("---------|----------- \n")
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   for p = 2:lenp
      @printf(" %1.1e | %4.2e  \n", H[p-1], errors[p-1])
   end
   @printf("---------|----------- \n")
   return dEh
   if minimum(errors) <= 1e-3 * maximum(errors)
      if verbose
         println("passed")
      end
      return true
   else
      @warn("""It seems the finite-difference test has failed, which indicates
      that there is an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
      return false
   end
end

function fdtest_vfun(at::Atoms{Float64}; verbose=true)
   X0 = positions(at)
   atX = mat(X0)[:,:]
   Xm = atX + 1e-2*rand(size(atX)...)
   X = vecs(Xm)[:]
   set_positions!(at, X)

   nlist = neighbourlist(at, cutoff(at))
   tmp = JuLIP.alloc_temp_d(at.calc, at) # tmp.R and tmp.dV
   _, R = neigs!(tmp.R, nlist, 1) # fillin tmp.R
   Rm = mat(R)[:,:]
   V = Potentials.evaluate!(tmp, at.calc, R)#/det(cell(at))/2
   Potentials.evaluate_d!(tmp.dV, tmp, at.calc, R)
   dV = -mat(tmp.dV)[:,:]

   errors = Float64[]
   # loop through finite-difference step-lengths
   @printf("---------|----------- \n")
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   Rhm = copy(Rm)
   for p = 2:11
      h = 0.1^p
      dVh = zero(dV)
      for n = 1:length(Rhm)
         Rhm[n] += h
         Rh = vecs(Rhm)[:]
         Vh = Potentials.evaluate!(tmp, at.calc, Rh)#/det(cell(at))/2
         dVh[n] = (V - Vh) / h
         Rhm[n] -= h
      end
      push!(errors, norm(dV - dVh, Inf))
      @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   @printf("---------|----------- \n")
   if minimum(errors) <= 1e-3 * maximum(errors)
      if verbose
         println("passed")
      end
      return true
   else
      @warn("""It seems the finite-difference test has failed, which indicates
      that there is an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
      return false
   end
end

function Wcb_vfun(at::Atoms{Float64})
   nlist = neighbourlist(at, cutoff(at))
	tmp = JuLIP.alloc_temp_d(at.calc, at) # tmp.R and tmp.dV
	_, R = neigs!(tmp.R, nlist, 1) # fillin tmp.R
	# Rref = mat(tmp.R)[:,:]
   V = Potentials.evaluate!(tmp, at.calc, R)/det(cell(at))/2

   # V = Potentials.evaluate!(tmp, at.calc, R)/det(cell(at))/2 - W0
	Potentials.evaluate_d!(tmp.dV, tmp, at.calc, R) # fillin tmp.dV 
   dV = -mat(tmp.dV)[:,:] 
   return V, dV
end

function fdtest(F::Function, dF::Function, at::Atoms{Float64}; verbose=true)
   H = Float64[]
   errors = Float64[]

   x0 = positions(at) |> mat
   x = deepcopy(x0)
   E = F(at)
   dE = dF(at)
   plen=4
   # loop through finite-difference step-lengths
   dEh = copy(dE)
   ldE = size(dE, 2)
   prg = Progress(3*(plen-1)*ldE)
   for p = 2:plen
      h = 0.1^p
      push!(H, h)
      
      for n = 1:size(dE,2)
           for j = 1:3
               x[j,n] += h
               set_positions!(at, x)
               Eh = F(at)
               # dEh[j,n] = (Eh - E) / h
               dEh[j,n] = (E - Eh)/h
               x[j,n] -= h
               set_positions!(at, x)
               ProgressMeter.update!(prg,3*(p-2)*ldE + 3*(n-1)+j)
           end
      end
      push!(errors, norm(dE - dEh, Inf))
   end
   @printf("---------|----------- \n")
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   for p = 2:plen
      @printf(" %1.1e | %4.2e  \n", H[p-1], errors[p-1])
   end
   @printf("---------|----------- \n")
   if minimum(errors) <= 1e-3 * maximum(errors)
      if verbose
         println("passed")
      end
      return true
   else
      @warn("""It seems the finite-difference test has failed, which indicates
      that there is an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
      return false
   end
end

function fdtest(F::Function, dF::Function, at::Atoms{Float64}, IB::Vector{Int64}; verbose=true)
   H = Float64[]
   errors = Float64[]

   x = positions(at) |> mat
   E = F(at)
   dE = dF(at)

   # loop through finite-difference step-lengths
   ldE = size(x, 2)
   prg = Progress(3*10*ldE)
   for p = 2:11
      h = 0.1^p
      push!(H, h)
      dEh = copy(dE)
      for n = 1:size(x,2)
         if !(n in IB)
            continue
         end 
         for j = 1:3
            x[j,n] += h
            set_positions!(at, x)
            Eh = F(at)
            dEh[j,n] = (Eh - E) / h
            x[j,n] -= h
            set_positions!(at, x)
            ProgressMeter.update!(prg,3*(p-2)*ldE + 3*(n-1)+j)
         end
      end
      push!(errors, norm(dE - dEh, Inf))
   end
   @printf("---------|----------- \n")
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   for p = 2:11
      @printf(" %1.1e | %4.2e  \n", H[p-1], errors[p-1])
   end
   @printf("---------|----------- \n")
   if minimum(errors) <= 1e-3 * maximum(errors)
      if verbose
         println("passed")
      end
      return true
   else
      @warn("""It seems the finite-difference test has failed, which indicates
      that there is an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
      return false
   end
end