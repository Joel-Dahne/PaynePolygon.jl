import PaynePolygon: stiffness_matrix

@testset "edges_to_triangles" begin
    edge_to_triangle1, edge_to_triangle2 = PaynePolygon.edges_to_triangles(36)
    # Compare to output from Matlab
    vars = matread("edges_to_triangles.mat")

    @test all(edge_to_triangle1 .== vars["EdgeToTriangle1"]')
    @test all(edge_to_triangle2 .== vars["EdgeToTriangle2"]')
end

@testset "boundary_and_interior" begin
    b, i, v = PaynePolygon.boundary_and_interior(9, 4, 3)
    @test b == BitSet([18, 19, 23, 24, 28, 29, 33, 34, 37, 38, 48, 49])
    @test i == BitSet([17, 25, 26, 27, 35, 36])
    @test v == BitSet([19, 23, 29, 33, 37, 49])
end

@testset "stiffness_matrix" begin
    # Check to that the returned matrix actually is hermitian
    @test ishermitian(stiffness_matrix(9, return_hermitian = false))
    @test ishermitian(stiffness_matrix(9, 4, 3, return_hermitian = false))

    # Check that it returns correct types
    @test stiffness_matrix(Float64, 9, return_hermitian = false) isa Matrix{Float64}
    @test stiffness_matrix(Float64, 9, return_hermitian = true) isa
          Hermitian{Float64,Matrix{Float64}}

    @test stiffness_matrix(Arb, 9, return_hermitian = false, return_arbmatrix = true) isa
          ArbMatrix
    @test stiffness_matrix(Arb, 9, return_hermitian = false, return_arbmatrix = false) isa
          Matrix{Arb}
    @test stiffness_matrix(Arb, 9, return_hermitian = true, return_arbmatrix = true) isa
          Hermitian{Arb,ArbMatrix}
    @test stiffness_matrix(Arb, 9, return_hermitian = true, return_arbmatrix = false) isa
          Hermitian{Arb,Matrix{Arb}}

    # Check that it runs when d and h are non-zero, doesn't check if
    # the result is correct
    @test PaynePolygon.stiffness_matrix(9, 4, 3) isa Hermitian{Float64}
    @test PaynePolygon.stiffness_matrix(Float64, 9, 4, 3) isa Hermitian{Float64}

    # Check that it gives reasonable results for N = 9, d = 4, h = 3
    # by comparing to precomputed results
    M = PaynePolygon.stiffness_matrix(9, 4, 3)
    @test eigvals(M, 1:5) ≈ [26.5765, 66.2828, 66.2828, 97.0068, 97.2862] rtol = 1e-5
end
