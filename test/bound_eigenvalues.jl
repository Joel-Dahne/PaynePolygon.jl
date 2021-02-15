@testset "separate_eigenvalues" begin
    M = Hermitian(Diagonal([1.0, 2.0, 3.0, 4.0]))

    @test PaynePolygon.separate_eigenvalues(M, 0) == -Inf
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4) == Inf

    Q = eigen(convert(Hermitian{Float64,Matrix{Float64}}, M)).vectors

    @test PaynePolygon.separate_eigenvalues(M, 0; Q) == -Inf
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1; Q) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2; Q) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3; Q) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4; Q) == Inf

    M = convert(Hermitian{Arb,Matrix{Arb}}, M)

    @test PaynePolygon.separate_eigenvalues(M, 0) == -Inf
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4) == Inf

    @test PaynePolygon.separate_eigenvalues(M, 0; Q) == -Inf
    @test 1 < PaynePolygon.separate_eigenvalues(M, 1; Q) <= 2
    @test 2 < PaynePolygon.separate_eigenvalues(M, 2; Q) <= 3
    @test 3 < PaynePolygon.separate_eigenvalues(M, 3; Q) <= 4
    @test PaynePolygon.separate_eigenvalues(M, 4; Q) == Inf

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

@testset "symtri!" begin
    M = PaynePolygon.stiffness_matrix(9, 4, 3)

    @test GenericLinearAlgebra.symtri!(copy(M)) == PaynePolygon.symtri!(copy(M))
end

@testset "_Array" begin
    M = PaynePolygon.stiffness_matrix(9, 4, 3)
    A = PaynePolygon.symtri!(M)

    @test Array(A.Q) == PaynePolygon._Array(A.Q)
end
