using PaynePolygon, Test, MAT, LinearAlgebra, Arblib, Random, GenericLinearAlgebra

@testset "PaynePolygon" begin
    include("stiffness_matrix.jl")
    include("bound_eigenvalues.jl")
end
