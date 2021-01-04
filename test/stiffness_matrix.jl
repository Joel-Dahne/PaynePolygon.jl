@testset "edges_to_triangles" begin
    edge_to_triangle1, edge_to_triangle2 = PaynePolygon.edges_to_triangles(36)
    # Compare to output from Matlab
    vars = matread("edges_to_triangles.mat")

    @test all(edge_to_triangle1 .== vars["EdgeToTriangle1"]')
    @test all(edge_to_triangle2 .== vars["EdgeToTriangle2"]')
end

@testset "stiffness_matrix" begin
    M = PaynePolygon.stiffness_matrix(return_hermitian = false)

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
end
