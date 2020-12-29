using PaynePolygon, Test, MAT, LinearAlgebra

@testset "PaynePolygon" begin
    include("stiffness_matrix.jl")
    include("bound_eigenvalues.jl")
end
