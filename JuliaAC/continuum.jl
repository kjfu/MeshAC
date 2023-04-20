"""
## module Continuum 

### Summary

This module contains the geometry info and energy density functional (Cauchy-Born) 
"""
module Continuum 

using LinearAlgebra
using JuLIP
# using JuLIP: AbstractAtoms, AbstractCalculator,
# 			 JVec, mat, vec, JMat, SVec, vecs, SMat

using JuLIP.Potentials #: evaluate_d!

import JuLIP: energy, forces, hessian_pos, hessian,
			  energy!, forces!


"""
`JData`: datatype for storing any data

some data which needs to be updated if the configuration (positions only!) has
changed too much.
"""
mutable struct JData{T<:AbstractFloat}
   max_change::T     # how much X may change before recomputing
   accum_change::T   # how much has it changed already
   data::Any
end

mutable struct Wcb{T <: AbstractFloat}
	W::Float64
	dW::Array{Float64,2}
	volT::Float64
end

abstract type AbstractContinuum{T} end

mutable struct Continuums{T} <: AbstractContinuum{T} #where T # <: AbstractFloat
	X::Vector{JVec{T}}			# nodes' positions
	nX::Int64					# number of nodes
	T::Vector{JVec{Int64}}		# elements
	nT::Int64					# number of elements
	# F::Vector{Vector{Int64}}	# face of elements
	# CN::Vector{Vector{Int64}}	# neighbours
	J::Vector{JMat{T}}			# Jacobian?
	DM::Vector{JMat{T}}			# deformation matrices
	## Wcb (energy density functional) related
	# calc::Wcb{T}	# calculator
	at::AbstractAtoms{T}			# unitCell with the same type of atomistic
	# C0::Vector{JMat{T}}			# original cell matrix
	V0::Vector{T}				# volumes of original cells
	data::Dict{Any,JData{T}}
	# E::T
	# Frc::Vector{JVec{T}}
	## Boundary conditions
	iFree::Vector{Int64}		# index of free nodes
end 

const I3 = Matrix(1.0*I,3,3)

function Wcb(idxT::Int64, Tet::Array{Int64,2}, X::Array{T,2}, U::Array{T,2}, at::Atoms{T}; V0=det(cell(at))/2) where {T}
	# at0 = deepcopy(at)
	t = Tet[1:4,idxT]
	# @show idxT
	# xx = mat(X[t])[:]
	# uu = mat(U[t])[:]
	xx = X[1:3,t]
	uu = U[1:3,t]
	
	J = zeros(T, 3, 3)
	Du = zeros(T, 3, 3)

	for j = 1:3
		J[:,j] = xx[:,j+1] - xx[:,1]
		Du[:,j] = uu[:,j+1] - uu[:,1]
	end

	∇u = (Du + J) / J # check Du or Du+J would be suitable for minimise!


	volT = abs(det(J)/6)

	nlist = neighbourlist(at, cutoff(at))
	tmp = JuLIP.alloc_temp_d(at.calc, at) # tmp.R and tmp.dV
	_, R = neigs!(tmp.R, nlist, 1) # fillin tmp.R
	Rref = mat(tmp.R)[:,:]
	W0 = Potentials.evaluate!(tmp, at.calc, R)/V0

	X0 = positions(at)
	C0 = cell(at)
	# @show idxT
	# @show cell(at)
	apply_defm!(at, ∇u)
	Rdef = ∇u * Rref
	Ridx = findall(x->x≤cutoff(at), vec(mapslices(x->norm(x,2), Rdef, dims=1)))
	# Ridx = findall(x->x<cutoff(at), vec(mapslices(x->norm(x,2), Rdef, dims=1)))
	Rref = Rref[:,Ridx]

	nlist = neighbourlist(at, cutoff(at))
	tmp = JuLIP.alloc_temp_d(at.calc, at) # tmp.R and tmp.dV
	_, R = neigs!(tmp.R, nlist, 1)

	W = Potentials.evaluate!(tmp, at.calc, R)/V0 - W0
	Potentials.evaluate_d!(tmp.dV, tmp, at.calc, R) # fillin tmp.dV 
	dV = mat(tmp.dV)[:,:]
	if size(dV,2) ≠ size(Rref,2)
		@show idxT
		@show uu
	end
	dW = dV*adjoint(Rref)/adjoint(J)
	# at = deepcopy(at0)
	set_positions!(at, X0)
	set_cell!(at, C0)
	return W, dW, volT
end # Wcb
# constructor

# function energy(contm::Continuum{T}) where {T}
function energy_frc(Tet::Array{Int64, 2}, Vex::Array{T,2}, U::Array{T,2}, at::Atoms{T}) where T
	E = zero(T)
	Frc = zeros(T, 3, size(Vex,2))
	for j = 1:size(Tet,2)
		# v0 = det(cell(contm.at))
		# E += energy(contm.calc, set_cell!(contm.at, (DM[j] * contm.C0[j]')') ) / v0 * contm.V0[j]
		W, dW, volT = Wcb(j, Tet, Vex, U, at)
		
		E += volT*W

		t = Tet[1:4,j]
		# Frc[:,t[1]] = -volT*sum(dW, dims=2)
		for n = 1:length(t)-1
			Frc[:,t[n+1]] += volT*dW[:,n]
			Frc[:,t[1]] -= volT*dW[:,n]
		end
	end
	Vbdry = findall(x->x ≠ 0, Vex[4,:])
	Frc[:,Vbdry] .= 0
	return E, Frc
end # energy_frc


end # module
