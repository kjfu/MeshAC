using JuLIPMaterials
using JuLIP

CB = JuLIPMaterials.CauchyBorn

at = bulk(:Fe)
r0 = rnn(:Fe)
calc = LennardJones(σ = r0) * C2Shift(2.7*r0)
set_calculator!(at, calc)
I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

W = CB.Wcb(at, calc, normalise = :atoms)
@show energy(at) - W(I)
@show CB.grad(W, I) - (-virial(at))

at = bulk(:Si)
set_calculator!(at, StillingerWeber())
I = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
W = CB.Wcb(at)
@show energy(at) - W(I)
@show CB.grad(W, I) - (-virial(at))
