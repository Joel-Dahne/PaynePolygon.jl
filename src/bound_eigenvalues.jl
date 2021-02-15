"""
    separate_eigenvalues(M::Hermitian, k::Integer; Q = nothing)

Compute a number `Λ` which separates the first `k` eigenvalues of `M`
from the rest.

`Q` should contain the (approximate) eigenvectors of `M` as columns,
or if `nothing` is given then it computes the eigenvectors by first
converting `M` to `Float64` and then computing the eigenvalues with
`eigen`. `Q` is then in both cases converted to a `Matrix{eltype(M)}`.

"""
function separate_eigenvalues(
    M::Hermitian,
    k::Integer;
    Q::Union{Nothing,AbstractMatrix} = nothing,
)
    N = size(M)[1]

    if k < 0
        throw(ArgumentError("k must be non-negative, not $k"))
    end
    if k == 0
        return convert(eltype(M), -Inf)
    end
    if k >= N
        return convert(eltype(M), Inf)
    end

    if isnothing(Q)
        # Compute approximate eigenvectors
        Q = eigen(convert(Hermitian{Float64,Matrix{Float64}}, M)).vectors
    else
        size(M) == size(Q) ||
            throw(ArgumentError("M and Q should have the same dimensions"))
    end

    # Convert Q to the appropriate type
    Q = convert(Matrix{eltype(M)}, Q)

    # Compute s from Lemma 2.4 in GS-O
    s = maximum(abs, Q' * Q - I)

    # Compute D
    MQ = M * Q
    D = Q' * MQ

    # Compute error bounds for D
    bounds = similar(D)

    Mv_norms = [norm(v) for v in eachcol(MQ)]
    M_norm = norm(M)
    for i = 1:N
        for j = 1:N
            bounds[i, j] = sqrt(3s) * (Mv_norms[i] + Mv_norms[j]) + 4s * M_norm
        end
    end

    # Compute the radii of the disks in the Gershgorin circle theorem
    Rs = similar(Mv_norms)
    for i in eachindex(Rs)
        Rs[i] = sum(abs(D[i, j]) + bounds[i, j] for j in [1:i-1; i+1:N])
    end

    # Compute upper bound of the first k disks
    Λ_lower = maximum([D[i, i] + Rs[i] for i = 1:k])

    # Compute lower bound of the remaining disks
    Λ_upper = minimum([D[i, i] - Rs[i] for i = k+1:N])

    if !(Λ_lower < Λ_upper)
        throw(
            ErrorException(
                "could not separate the first $k eigenvalues from the rest," *
                " got lower bound $Λ_lower and upper bound $Λ_upper",
            ),
        )
    end

    # We want the bound to be as large as possible, so take the upper bound
    return Λ_upper
end

"""
    separate_eigenvalues(M::Hermitian{Arb}, k::Integer; Q = nothing)

Compute a number `Λ` which separates the first `k` eigenvalues of `M`
from the rest.

The number is computed in a rigorous way as to guarantee the
separation.

`Q` should contain the (approximate) eigenvectors of `M` as columns,
or if `nothing` is given then it computes the eigenvectors by first
converting `M` to `Float64` and then computing the eigenvalues with
`eigen`. `Q` is then in both cases converted to an `ArbMatrix`.

"""
function separate_eigenvalues(
    M::Hermitian{Arb},
    k::Integer;
    Q::Union{Nothing,AbstractMatrix} = nothing,
)
    N = size(M)[1]

    if k < 0
        throw(ArgumentError("k must be non-negative, not $k"))
    end
    if k == 0
        return convert(eltype(M), -Inf)
    end
    if k >= N
        return convert(eltype(M), Inf)
    end

    if isnothing(Q)
        # Compute approximate eigenvectors
        Q = eigen(convert(Hermitian{Float64,Matrix{Float64}}, M)).vectors
    else
        size(M) == size(Q) ||
            throw(ArgumentError("M and Q should have the same dimensions"))
    end

    # Convert Q to an ArbMatrix if needed
    if !(Q isa ArbMatrix)
        Q = ArbMatrix(Q)
    end
    Q_transpose = Arblib.transpose!(similar(Q), Q)

    # M as an ArbMatrix
    M = ArbMatrix(M)

    # Compute s from Lemma 2.4 in GS-O
    s = let ss = Q_transpose * Q
        s = zero(ss[1])
        for i in CartesianIndices(ss)
            if i.I[1] == i.I[2]
                s = Arblib.max!(s, s, abs(ss[i] - 1))
            else
                s = Arblib.max!(s, s, abs(ss[i]))
            end
        end
        s
    end

    # Compute D
    MQ = M * Q
    D = Q_transpose * MQ

    # Compute error bounds for D
    bounds = similar(D)

    # TODO: Check that norm produces rigorous results
    Mv_norms = [norm(v) for v in eachcol(MQ)]
    M_norm = Arblib.frobenius_norm!(Arb(), M)

    for i = 1:N
        for j = 1:N
            bounds[i, j] = sqrt(3s) * (Mv_norms[i] + Mv_norms[j]) + 4s * M_norm
        end
    end

    # Compute the radii of the disks in the Gershgorin circle theorem
    Rs = similar(Mv_norms)
    for i in eachindex(Rs)
        Rs[i] = sum(abs(D[i, j]) + bounds[i, j] for j in [1:i-1; i+1:N])
    end

    # Compute upper bound of the first k disks
    Λ_lower = maximum([D[i, i] + Rs[i] for i = 1:k])

    # Compute lower bound of the remaining disks
    Λ_upper = minimum([D[i, i] - Rs[i] for i = k+1:N])

    if !(Λ_lower < Λ_upper)
        throw(
            ErrorException(
                "could not separate the first $k eigenvalues from the rest," *
                " got lower bound $Λ_lower and upper bound $Λ_upper",
            ),
        )
    end

    # We want the bound to be as large as possible, so take the upper bound
    return Λ_upper
end

"""
    symtri!(A::Hermitian)

Reduce `A` to symmetric tridiagonal form. The code is the same as in
[GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/)
but with debug info and parallelized.

It has the following license

The MIT License (MIT)

Copyright (c) 2014-2018 Andreas Noack

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
function symtri!(A::Hermitian{T}) where {T<:Real}
    A.uplo == 'U' || throw(ArgumentError("only A.uplo = 'U' is supported, got $A.uplo"))
    return symtriUpper!(A.data)
end

# Assume that upper triangle stores the relevant part
function symtriUpper!(
    AS::StridedMatrix{T},
    τ = zeros(T, size(AS, 1) - 1),
    u = Vector{T}(undef, size(AS, 1)),
) where {T}
    n = LinearAlgebra.checksquare(AS)

    # We ignore any non-real components of the diagonal
    @inbounds for i = 1:n
        AS[i, i] = real(AS[i, i])
    end

    @inbounds for k = 1:(n-2+!(T <: Real))
        # This is a bit convoluted method to get the conjugated vector
        # but conjugation is required for correctness of arrays of
        # quaternions. Eventually, it should be sufficient to write
        # vec(x') but it currently (July 10, 2018) hits a bug in
        # LinearAlgebra
        τk = LinearAlgebra.reflector!(vec(transpose(view(AS, k, (k+1):n)')))
        τ[k] = τk'

        Threads.@threads for j = (k+1):n
            tmp = AS[k+1, j]
            for i = (k+2):j
                tmp += AS[k, i] * AS[i, j]
            end
            for i = (j+1):n
                tmp += AS[k, i] * AS[j, i]'
            end
            u[j] = τk' * tmp
        end

        vcAv = u[k+1]
        for i = (k+2):n
            vcAv += u[i] * AS[k, i]'
        end
        ξ = real(vcAv * τk)

        Threads.@threads for j = (k+1):n
            ujt = u[j]
            hjt = j > (k + 1) ? AS[k, j] : one(ujt)
            ξhjt = ξ * hjt
            for i = (k+1):(j-1)
                hit = i > (k + 1) ? AS[k, i] : one(ujt)
                AS[i, j] -= hit' * ujt + u[i]' * hjt - hit' * ξhjt
            end
            AS[j, j] -= 2 * real(hjt' * ujt) - abs2(hjt) * ξ
        end
    end
    GenericLinearAlgebra.SymmetricTridiagonalFactorization(
        'U',
        AS,
        τ,
        SymTridiagonal(real(diag(AS)), real(diag(AS, 1))),
    )
end

"""
    _Array!(Q::GenericLinearAlgebra.EigenQ)

Compute `Array(Q)` in a more efficient way. The code is similar to
[GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/)
but optimized and parallelized.

It has the following license

The MIT License (MIT)

Copyright (c) 2014-2018 Andreas Noack

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
function _Array(Q::GenericLinearAlgebra.EigenQ)
    Q.uplo == 'U' || throw(ArgumentError("only Q.uplo = 'U' supported."))

    n = size(Q, 1)
    B = Matrix{eltype(Q)}(I, n, n)

    for k = length(Q.τ):-1:1
        for j = 1:n
            b = view(B, :, j)
            vb = B[k+1, j]
            for i = (k+2):n
                vb += Q.factors[k, i] * B[i, j]
            end
            τkvb = Q.τ[k]' * vb
            B[k+1, j] -= τkvb
            for i = (k+2):n
                B[i, j] -= Q.factors[k, i]' * τkvb
            end
        end
    end

    return B
end
