using FEMBasis
using JuLIP


# tri-linear basis
code = FEMBasis.create_basis(
    :Cube8,
    "8 node tri-linear cubic element produced by DealII",
    (
        (0.0, 0.0, 0.0), # N1
        (1.0, 0.0, 0.0), # N2
        (1.0, 0.0, 1.0), # N3
        (0.0, 0.0, 1.0), # N4
        (0.0, 1.0, 0.0), # N5
        (1.0, 1.0, 0.0), # N6
        (1.0, 1.0, 1.0), # N7
        (0.0, 1.0, 1.0), # N8
    ),
    :(1 + u + v + w + u*v + v*w + w*u + u*v*w),
)

eval(code)

# codeh = FEMBasis.create_basis(
#     :Cubeh8,
#     "8 node tri-linear cubic element produced by DealII",
#     (
#         (0.0,  0.0,  0.0), # N1
#         (0.25, 0.0,  0.0), # N2
#         (0.25, 0.0,  0.25), # N3
#         (0.0,  0.0,  0.25), # N4
#         (0.0,  0.25, 0.0), # N5
#         (0.25, 0.25, 0.0), # N6
#         (0.25, 0.25, 0.25), # N7
#         (0.0,  0.25, 0.25), # N8
#     ),
#     :(1 + u + v + w + u*v + v*w + w*u + u*v*w),
# )
# eval(codeh)

ξ = zeros(Float64, 8, 1)
x = rand(3,1)
h = 0.25
eval_basis!(Cube8, ξ, x/h)

dξ = zeros(Float64, 3, 1)
eval_dbasis!(Cube8, dξ, x/h)
dξ = dξ/h