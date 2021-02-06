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
    ss = abs.(Q' * Q - I)
    s = maximum(ss)

    # Compute D
    MQ = M * Q
    D = Q' * MQ

    # Compute error bounds for D
    bounds = similar(D)

    Mv_norms = [norm(v) for v in eachcol(MQ)]
    M_norm = norm(M)

    for (i, vi) in enumerate(eachcol(Q))
        for (j, vj) in enumerate(eachcol(Q))
            bounds[i, j] = sqrt(3s) * (Mv_norms[i] + Mv_norms[j]) + 4s * M_norm
        end
    end

    # Compute the radii of the disks in the Gershgorin circle theorem
    Rs = similar(Mv_norms)
    for i in eachindex(Rs)
        Rs[i] = sum(D[i, j] + bounds[i, k] for j in [1:i-1; i+1:N])
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
    ss = abs.(Q_transpose * Q - I)
    s = maximum(ss)

    # Compute D
    MQ = M * Q
    D = Q_transpose * MQ

    # Compute error bounds for D
    bounds = similar(D)

    # TODO: Check that norm produces rigorous results
    Mv_norms = [norm(v) for v in eachcol(MQ)]
    M_norm = Arblib.frobenius_norm!(Arb(), M)

    for (i, vi) in enumerate(eachcol(Q))
        for (j, vj) in enumerate(eachcol(Q))
            bounds[i, j] = sqrt(3s) * (Mv_norms[i] + Mv_norms[j]) + 4s * M_norm
        end
    end

    # Compute the radii of the disks in the Gershgorin circle theorem
    Rs = similar(Mv_norms)
    for i in eachindex(Rs)
        Rs[i] = sum(D[i, j] + bounds[i, k] for j in [1:i-1; i+1:N])
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
