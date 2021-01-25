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
    M = PaynePolygon.stiffness_matrix(36, return_hermitian = false)

    # Check that it's Hermitian
    @test ishermitian(M)

    # Compare with some values computed from the Matlab code
    @test sum(M) == 2239488
    @test sum(!iszero, M) == 57348

    M = PaynePolygon.stiffness_matrix(18, return_hermitian = true)
    @test M isa Hermitian

    # Check that it returns correct types
    @test eltype(PaynePolygon.stiffness_matrix(9)) == Int
    @test eltype(PaynePolygon.stiffness_matrix(Float64, 9)) == Float64
    @test eltype(PaynePolygon.stiffness_matrix(Arb, 9)) == Arb
    @test eltype(PaynePolygon.stiffness_matrix(9, return_hermitian = false)) == Int
    @test eltype(PaynePolygon.stiffness_matrix(Float64, 9, return_hermitian = false)) ==
          Float64
    @test eltype(PaynePolygon.stiffness_matrix(Arb, 9, return_hermitian = false)) == Arb

    # Check that it runs when d and h are non-zero, doesn't check if
    # the result is correct
    @test PaynePolygon.stiffness_matrix(9, 4, 3) isa Hermitian{Int}
    @test PaynePolygon.stiffness_matrix(Float64, 9, 4, 3) isa Hermitian{Float64}
end
