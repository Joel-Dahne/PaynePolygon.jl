module PaynePolygon

using LinearAlgebra,
    Arblib,
    ArbTools,
    JLD,
    MethodOfParticularSolutions,
    Nemo,
    Plots,
    GenericLinearAlgebra,
    Printf

include("stiffness_matrix.jl")
include("bound_eigenvalues.jl")

include("plotting.jl")

include("handle_eigenfunctions.jl")
include("bounded_by.jl")

end # module
