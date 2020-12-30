@testset "separate_eigenvalues" begin
    M = Hermitian(Diagonal([1.0, 2.0, 3.0, 4.0]))

    @test PaynePolygon.separate_eigenvalues(M, 0) == -Inf
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4) == Inf

    M = Hermitian(rand(Random.MersenneTwister(42), 4, 4))
    M_arb = convert(Hermitian{Arb,Matrix{Arb}}, M)
    @test PaynePolygon.separate_eigenvalues(M_arb, 0) == -Inf
    @test abs(
        PaynePolygon.separate_eigenvalues(M, 1) -
        PaynePolygon.separate_eigenvalues(M_arb, 1),
    ) <= 1e-6
    @test abs(
        PaynePolygon.separate_eigenvalues(M, 2) -
        PaynePolygon.separate_eigenvalues(M_arb, 2),
    ) <= 1e-6
    @test abs(
        PaynePolygon.separate_eigenvalues(M, 3) -
        PaynePolygon.separate_eigenvalues(M_arb, 3),
    ) <= 1e-6
    @test PaynePolygon.separate_eigenvalues(M, 4) == Inf
end
