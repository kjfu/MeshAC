
# JuLIPAC.jl master file

module JuLIPAC

# the atoms are build using JuLIP
using JuLIP

include("continuum.jl")

include("AtC.jl")

include("dihole.jl")

include("Solve.jl")

include("adaptive.jl")

include("FIO.jl")

include("Testing.jl")

include("utils.jl")

end # module
