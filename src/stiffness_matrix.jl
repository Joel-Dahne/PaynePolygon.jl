"""
    edges_to_triangle(N::Integer = 9)

Compute edges between triangles for a uniform triangular grid with
side lengths `1/N`. Used in [`stiffness_matrix`](@ref). The result is
exact.
"""
function edges_to_triangles(N::Integer = 9)
    function put!(v1, v2, x, c)
        if x > length(v1)
            x -= length(v1)
        end
        if x < 0
            x += length(v1)
        end

        if v1[x] == 0
            v1[x] = c
        elseif v2[x] == 0
            v2[x] = c
        else
            throw(ErrorException("something strange going on"))
        end
    end

    num_edges = div(3N^2 - N, 2)

    edge_to_t1 = zeros(Int, 6num_edges)
    edge_to_t2 = zeros(Int, 6num_edges)

    current_triangle = 1

    for i = 1:6
        for row = 1:N-1
            for index = 1:2row-2
                if index % 2 == 1
                    x =
                        (i - 1) * num_edges +
                        div(3(row - 1)^2 + (row - 1), 2) +
                        3div(index - 1, 2)
                    put!(edge_to_t1, edge_to_t2, x + 1, current_triangle)
                    put!(edge_to_t1, edge_to_t2, x + 2, current_triangle)
                    put!(edge_to_t1, edge_to_t2, x + 3, current_triangle)
                else
                    x =
                        (i - 1) * num_edges +
                        div(3(row - 1)^2 + (row - 1), 2) +
                        3(div(index, 2) - 1)
                    put!(edge_to_t1, edge_to_t2, x + 3, current_triangle)
                    put!(edge_to_t1, edge_to_t2, x + 4, current_triangle)
                    put!(
                        edge_to_t1,
                        edge_to_t2,
                        (i - 1) * num_edges +
                        div(3(row - 2)^2 + (row - 2), 2) +
                        3div(index - 2, 2) +
                        2,
                        current_triangle,
                    )
                end
                current_triangle += 1
            end

            # index = 2row - 1
            x = (i - 1) * num_edges + div(3 * row^2 + row, 2)
            put!(
                edge_to_t1,
                edge_to_t2,
                (i - 1) * num_edges + div(3 * row^2 + row, 2) - 1,
                current_triangle,
            )
            put!(
                edge_to_t1,
                edge_to_t2,
                (i - 1) * num_edges + div(3 * row^2 + row, 2),
                current_triangle,
            )
            put!(
                edge_to_t1,
                edge_to_t2,
                (i - 1) * num_edges +
                div(3N^2 - N, 2) +
                div(3(row - 1)^2 + (row - 1), 2) +
                1,
                current_triangle,
            )

            current_triangle += 1
        end

        for index = 1:2N-2
            if index % 2 == 1
                x = (i - 1) * num_edges + div(3(N - 1)^2 + (N - 1), 2) + index
                put!(edge_to_t1, edge_to_t2, x, current_triangle)
                put!(edge_to_t1, edge_to_t2, x + 1, current_triangle)
            else
                x = (i - 1) * num_edges + div(3(N - 1)^2 + (N - 1), 2)
                put!(edge_to_t1, edge_to_t2, x + index, current_triangle)
                put!(edge_to_t1, edge_to_t2, x + index + 1, current_triangle)
                put!(
                    edge_to_t1,
                    edge_to_t2,
                    (i - 1) * num_edges +
                    div(3(N - 2)^2 + (N - 2), 2) +
                    3div(index - 2, 2) +
                    2,
                    current_triangle,
                )
            end
            current_triangle += 1
        end

        x = (i - 1) * num_edges + div(3N^2 - N, 2)
        put!(edge_to_t1, edge_to_t2, x, current_triangle)
        put!(edge_to_t1, edge_to_t2, x + div(3(N - 1)^2 + (N - 1), 2) + 1, current_triangle)
        current_triangle += 1
    end

    return edge_to_t1, edge_to_t2
end

""""
    boundary_and_interior(N::Integer = 9, d::Integer = 4, h::Integer = 3)

Return three bitsets, the first containing all triangles intersecting
the boundary of the interior domain, the second containing all
triangles in their interior and the third containing all triangles
which intersect the vertices of the interior.

Only returns the indices in the fundamental domain.
"""
function boundary_and_interior(N::Integer = 9, d::Integer = 4, h::Integer = 3)
    h % 3 == 0 || throw(ArgumentError("h must be divisible by 3, got h = $h"))

    boundary = BitSet()
    interior = BitSet()
    vertices = BitSet()

    d == 0 && h == 0 && return boundary, interior, vertices

    # The height of the interior domain is h, iterate over each row in
    # the height.
    for i = 1:h
        # Determine the index for the first triangle on the row and
        # the number of triangles on the row before we hit the
        # boundary of the interior domain.
        if i <= h ÷ 3
            row_interior_start = (d + i - 1)^2 + 1
            row_interior_length = 4(i - 1) + 1
        else
            row_interior_start = (d + i - 1)^2 + 1
            row_interior_length = 2(h - i)
        end

        # Add the indices for the triangles in the interior
        for j = 1:row_interior_length
            push!(interior, row_interior_start + j - 1)

            push!(interior, row_interior_start + 2(d + i) - j - 1)
        end

        # Add the indices for the triangles on the boundary
        push!(boundary, row_interior_start + row_interior_length)
        push!(boundary, row_interior_start + row_interior_length + 1)

        push!(boundary, row_interior_start + 2(d + i - 1) - row_interior_length)
        push!(boundary, row_interior_start + 2(d + i - 1) - row_interior_length - 1)
    end

    # 60 degree vertices
    row = d + h ÷ 3 # Which row the vertex is on
    column = 4h ÷ 3 - 1 # How far in to move

    push!(vertices, (row - 1)^2 + column)
    push!(vertices, row^2 + column + 1)

    push!(vertices, row^2 - column + 1)
    push!(vertices, (row + 1)^2 - column)

    # 30 degree vertices
    push!(vertices, (d + h - 1)^2 + 1)
    push!(vertices, (d + h)^2)

    return boundary, interior, vertices
end


"""
    stiffness_matrix(N::Integer = 9, d::Integer = 0, h::Integer = 0; return_hermitian = true)
    stiffness_matrix(T, N::Integer = 9, d::Integer = 0, h::Integer = 0; return_hermitian = true)

Return the stiffness matrix for the triangulation of the problem.

If `T` is given convert the result to this type. The elements are
rational numbers so the result will be exact for `T<:Rational` if it
can hold sufficiently large integers (around `N^2`). If `T = Arb` then
the result will be an enclosure of the true result.

The parameter `N` determines the fines of the grid. `d` and `h`
determines the placement of the interior domains, setting them to zero
gives no interior domains.

The result is always hermitian, if `return_hermitian` is true then
return an explicitly hermitian matrix (of type `Hermitian`), otherwise
return a normal matrix.
"""
function stiffness_matrix(
    T::Type{<:Real},
    N::Integer = 9,
    d::Integer = 0,
    h::Integer = 0;
    return_hermitian = true,
)
    num_triangles = N^2 # Number of triangles in each fundamental domain
    edges = collect(zip(edges_to_triangles(N)...))
    boundary, interior, vertices = boundary_and_interior(N, d, h)

    # Remove edges we don't want
    # 1. Edges for which one triangle is in the interior of the interior domains
    # 2. Edges for triangles which meet at a vertex of the interior domains
    filter!(
        ((t1, t2),) ->
            !((t1 % num_triangles ∈ interior) || (t2 % num_triangles ∈ interior)) &&
                !((t1 % num_triangles ∈ vertices) && (t2 % num_triangles ∈ vertices)),
        edges,
    )

    num_edges = length(edges)
    stiffness_matrix = zeros(T, num_edges, num_edges)
    mass_diagonal = zeros(T, num_edges)

    # Notice that the stiffness matrix and the mass matrix both
    # contain the factor 1/sqrt(3). However in the end they all cancel
    # out so we remove them from the beginning.
    for i in eachindex(edges)
        for j in eachindex(edges)
            (T11, T21) = edges[i]
            (T12, T22) = edges[j]

            if min(T11, T21) != min(T12, T22) &&
               max(T11, T21) != max(T12, T22) &&
               min(T11, T21) != max(T12, T22) &&
               max(T11, T21) != min(T12, T22)
                # 0 hits, do nothing
                stiffness_matrix[i, j] = 0
            elseif min(T11, T21) == min(T12, T22) && max(T11, T21) == max(T12, T22)
                # 2 hits
                if i != j
                    throw(
                        ErrorException("Something strange going on at (i, j) = $((i, j))"),
                    )
                end
                if T11 % num_triangles ∈ boundary && T21 % num_triangles ∈ boundary
                    # Both of the triangles are on the boundary (1 - 1)
                    stiffness_matrix[i, j] = 4
                    mass_diagonal[i] = 1 // 4N^2
                elseif T11 % num_triangles ∈ boundary || T21 % num_triangles ∈ boundary
                    # One of the triangles are on the boundary (2 - 2)
                    stiffness_matrix[i, j] = 12
                    mass_diagonal[i] = 3 // 8N^2
                else
                    # Neither triangle is on the boundary (0 - 0)
                    stiffness_matrix[i, j] = 8
                    mass_diagonal[i] = 1 // 2N^2
                end
            else
                # 1 hit
                if (T11 % num_triangles ∈ boundary || T21 % num_triangles ∈ boundary) &&
                   (T12 % num_triangles ∈ boundary || T22 % num_triangles ∈ boundary)
                    # Both of the edges have one triangle on the boundary (1 - 2)
                    stiffness_matrix[i, j] = -2
                else
                    # (0 - 0) or (2 - 0)
                    stiffness_matrix[i, j] = -2
                end
            end
        end
    end

    mass_matrix = Diagonal(mass_diagonal)

    if return_hermitian
        return Hermitian(inv(mass_matrix) * stiffness_matrix)
    else
        return inv(mass_matrix) * stiffness_matrix
    end
end

stiffness_matrix(N::Integer = 9, d::Integer = 0, h::Integer = 0; return_hermitian = true) =
    stiffness_matrix(Float64, N, d, h; return_hermitian)
