function edges_to_triangles(N::Integer = 36)
    function put!(v1, v2, x, c, debug = "")
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
            throw(ErrorException("something strange going on at $debug"))
        end
    end

    num_edges = div(3N^2 - N, 2)
    num_triangles = (N - 1)^2

    edge_to_triangle1 = zeros(Int, 6num_edges)
    edge_to_triangle2 = zeros(Int, 6num_edges)

    current_triangle = 1

    for i = 1:6
        for row = 1:N-1
            for index = 1:2row-2
                if index % 2 == 1
                    x =
                        (i - 1) * num_edges +
                        div(3(row - 1)^2 + (row - 1), 2) +
                        3div(index - 1, 2)
                    put!(edge_to_triangle1, edge_to_triangle2, x + 1, current_triangle, "1")
                    put!(edge_to_triangle1, edge_to_triangle2, x + 2, current_triangle, "1")
                    put!(edge_to_triangle1, edge_to_triangle2, x + 3, current_triangle, "1")
                else
                    x =
                        (i - 1) * num_edges +
                        div(3(row - 1)^2 + (row - 1), 2) +
                        3(div(index, 2) - 1)
                    put!(edge_to_triangle1, edge_to_triangle2, x + 3, current_triangle, "2")
                    put!(edge_to_triangle1, edge_to_triangle2, x + 4, current_triangle, "2")
                    put!(
                        edge_to_triangle1,
                        edge_to_triangle2,
                        (i - 1) * num_edges +
                        div(3(row - 2)^2 + (row - 2), 2) +
                        3div(index - 2, 2) +
                        2,
                        current_triangle,
                        "2",
                    )
                end
                current_triangle += 1
            end

            # index = 2row - 1
            x = (i - 1) * num_edges + div(3 * (row)^2 + (row), 2)
            put!(
                edge_to_triangle1,
                edge_to_triangle2,
                (i - 1) * num_edges + div(3 * (row)^2 + (row), 2) - 1,
                current_triangle,
                "3",
            )
            put!(
                edge_to_triangle1,
                edge_to_triangle2,
                (i - 1) * num_edges + div(3 * (row)^2 + (row), 2),
                current_triangle,
                "3",
            )
            put!(
                edge_to_triangle1,
                edge_to_triangle2,
                (i - 1) * num_edges +
                div(3N^2 - N, 2) +
                div(3(row - 1)^2 + (row - 1), 2) +
                1,
                current_triangle,
                (N, row),
            )

            current_triangle += 1
        end

        for index = 1:2N-2
            if index % 2 == 1
                x = (i - 1) * num_edges + div(3(N - 1)^2 + (N - 1), 2) + index
                put!(edge_to_triangle1, edge_to_triangle2, x, current_triangle, "4")
                put!(edge_to_triangle1, edge_to_triangle2, x + 1, current_triangle, "4")
            else
                x = (i - 1) * num_edges + div(3(N - 1)^2 + (N - 1), 2)
                put!(edge_to_triangle1, edge_to_triangle2, x + index, current_triangle, "5")
                put!(
                    edge_to_triangle1,
                    edge_to_triangle2,
                    x + index + 1,
                    current_triangle,
                    "5",
                )
                put!(
                    edge_to_triangle1,
                    edge_to_triangle2,
                    (i - 1) * num_edges +
                    div(3(N - 2)^2 + (N - 2), 2) +
                    3div(index - 2, 2) +
                    2,
                    current_triangle,
                    (i, num_edges, N, index),
                )
            end
            current_triangle += 1
        end

        x = (i - 1) * num_edges + div(3N^2 - N, 2)
        put!(edge_to_triangle1, edge_to_triangle2, x, current_triangle, "6")
        put!(
            edge_to_triangle1,
            edge_to_triangle2,
            x + div(3(N - 1)^2 + (N - 1), 2) + 1,
            current_triangle,
            "6",
        )
        current_triangle += 1
    end

    return edge_to_triangle1, edge_to_triangle2
end

function stiffness_matrix(N::Integer = 36; return_hermitian = true)
    num_edges = div(3N^2 - N, 2)
    edge_to_triangle1, edge_to_triangle2 = edges_to_triangles(N)

    stiffness_matrix = zeros(6num_edges, 6num_edges)

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
                stiffness_matrix[i, j] = 8 / 3 * sqrt(3)
            else
                # 1 hit
                stiffness_matrix[i, j] = -2 / 3 * sqrt(3)
            end
        end
    end

    stiffness_matrix *= 3N^2 / 2 / (sqrt(3) / 4)

    if return_hermitian
        return Hermitian(stiffness_matrix)
    else
        return stiffness_matrix
    end
end
