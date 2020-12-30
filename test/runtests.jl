using PaynePolygon, Test, MAT, LinearAlgebra, Arblib, Random

@testset "PaynePolygon" begin
    include("stiffness_matrix.jl")
    include("bound_eigenvalues.jl")
end
