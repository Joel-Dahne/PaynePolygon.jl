"""
    edges_to_triangle(N::Integer = 36)

Compute edges between triangles, used in [`stiffness_matrix`](@ref).
All elements are integers so the result is exact.
"""
function edges_to_triangles(N::Integer = 36)
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
    boundary_and_interior(N::Integer, d::Integer, h::Integer)

Return two bitsets, the first containing all triangles intersecting
the boundary of the interior domain and the second containing all
triangles in their interior.

Only returns the indices in the fundamental domain.
"""
function boundary_and_interior(N::Integer, d::Integer, h::Integer)
    @assert mod(h, 3) == 0

    boundary = BitSet()
    interior = BitSet()

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

    return boundary, interior
end


"""
    stiffness_matrix(N::Integer = 36; return_hermitian = true)
    stiffness_matrix(T, N::Integer = 36; return_hermitian = true)

Return the stiffness matrix for the triangulation of the problem. All
the elements are integers so the result is exact.

If `T` is given convert the result to this type, the result will be
exact as long as `T` can represent `16N^2` and `-4N^2` exactly.

The parameter `N` determines the fines of the grid. The result is
always hermitian, if `return_hermitian` is true then return an
explicitly hermitian matrix (of type `Hermitian`), otherwise return a
normal matrix.
"""
function stiffness_matrix(T, N::Integer = 36; return_hermitian = true)
    num_edges = div(3N^2 - N, 2)
    edge_to_triangle1, edge_to_triangle2 = edges_to_triangles(N)

    stiffness_matrix = zeros(T, 6num_edges, 6num_edges)

    for i = 1:6num_edges
        for j = 1:6num_edges
            T11 = edge_to_triangle1[i]
            T12 = edge_to_triangle1[j]
            T21 = edge_to_triangle2[i]
            T22 = edge_to_triangle2[j]

            if min(T11, T22) != min(T12, T22) &&
               max(T11, T21) != max(T12, T22) &&
               min(T11, T21) != max(T12, T22) &&
               max(T11, T21) != min(T12, T22)
                # 0 hits, do nothing
                stiffness_matrix[i, j] = 0
            elseif min(T11, T21) == min(T12, T22) && max(T11, T21) == max(T12, T22)
                # 2 hits
                if i != j
                    throw(ErrorException("Something strange going on at (i, j) = $((i, j))"))
                end
                stiffness_matrix[i, j] = 4
            else
                # 1 hit
                stiffness_matrix[i, j] = -1
            end
        end
    end

    stiffness_matrix *= 4N^2

    if return_hermitian
        return Hermitian(stiffness_matrix)
    else
        return stiffness_matrix
    end
end

stiffness_matrix(N::Integer = 36; return_hermitian = true) =
    stiffness_matrix(Int, N; return_hermitian)
