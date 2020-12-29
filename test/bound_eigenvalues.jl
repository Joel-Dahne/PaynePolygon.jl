@testset "separate_eigenvalues" begin
    M = Hermitian(Diagonal([1.0, 2.0, 3.0, 4.0]))

    @test PaynePolygon.separate_eigenvalues(M, 0) == 0.0
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4) == Inf
end
