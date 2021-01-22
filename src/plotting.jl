"""
    mesh_triangle_midpoints(N::Integer = 36)

Return the midpoints of the triangles in the mesh. This is NOT
rigorous.

"""
function mesh_triangle_midpoints(N::Integer = 36)
    num_triangles = N^2

    xs = zeros(6num_triangles)
    ys = zeros(6num_triangles)

    # Height and side length of each triangle
    height = sqrt(3) / 2N
    side = 1 / N

    # Compute the x and y coordinates of the mesh-triangles in the
    # first big triangle.

    # x and y increment along each row
    x_inc = cospi(2 // 3) * side / 2
    y_inc = sinpi(2 // 3) * side / 2

    for row = 1:N
        row_start = (row - 1)^2

        # First part of the triangles

        # x and y coordinates for the first triangle in the row
        x₀ = side * (row - 1) + side / 2
        y₀ = side * tan(π / 6) / 2

        for i = 1:2:(2row-1)
            xs[row_start+i] = x₀ + (i - 1) * x_inc
            ys[row_start+i] = y₀ + (i - 1) * y_inc
        end

        # Second part of the triangles

        # x and y coordinates for the second triangle in the row
        x₀ = side * (row - 1)
        y₀ = height - side * tan(π / 6) / 2

        for i = 2:2:(2row-1)
            xs[row_start+i] = x₀ + (i - 2) * x_inc
            ys[row_start+i] = y₀ + (i - 2) * y_inc
        end
    end

    # Compute the x and y coordinates of the mesh-triangles in the
    # other triangles by rotation
    for triangle = 2:6
        θ = (triangle - 1) * π / 3
        s, c = sincos(θ)

        for i = 1:num_triangles
            x = xs[i]
            y = ys[i]
            xs[num_triangles*(triangle-1)+i] = x * c - y * s
            ys[num_triangles*(triangle-1)+i] = x * s + y * c
        end
    end

    return xs, ys
end

"""
    triangle_side(c1, c2)

Given `c1` and `c2` compute endpoints of the edge going between
triangles centered at these points.
"""
function triangle_edge(c1, c2)
    v = c2 - c1

    M = [cospi(1 // 3) -sinpi(1 // 3); sinpi(1 // 3) cospi(1 // 3)]

    # c1 + (v rotated π/3 positive direction)
    e1 = c1 + M * v
    # c1 + (v rotated -π/3 positive direction)
    e2 = c1 + inv(M) * v

    return e1, e2
end

function plot_mesh(N::Integer; excluded = BitSet([7]), marked = BitSet([8]))
    xs, ys = mesh_triangle_midpoints(N)
    edge_to_triangle1, edge_to_triangle2 = edges_to_triangles(N)

    num_triangles = N^2
    edges = collect(zip(edge_to_triangle1, edge_to_triangle2))
    indices = findall(
        ((t1, t2),) ->
            !((t1 % num_triangles ∈ excluded) || (t2 % num_triangles ∈ excluded)),
        edges,
    )
    edge_to_triangle1 = edge_to_triangle1[indices]
    edge_to_triangle2 = edge_to_triangle2[indices]

    xs_edges = zeros(2, length(edge_to_triangle1))
    ys_edges = zeros(2, length(edge_to_triangle1))
    colors = permutedims([
        ifelse(
            edge_to_triangle1[i] % num_triangles ∈ marked ||
            edge_to_triangle2[i] % num_triangles ∈ marked,
            :red,
            :blue,
        ) for i in eachindex(edge_to_triangle1)
    ])

    for (i, (t1, t2)) in enumerate(zip(edge_to_triangle1, edge_to_triangle2))
        e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]])
        xs_edges[1, i] = e1[1]
        xs_edges[2, i] = e2[1]
        ys_edges[1, i] = e1[2]
        ys_edges[2, i] = e2[2]
    end

    h = sqrt(3) / 2
    pl = plot(
        [1, 0.5, -0.5, -1, -0.5, 0.5, 1],
        [0, h, h, 0, -h, -h, 0],
        aspect_ratio = true,
        legend = :none,
    )

    scatter!(pl, xs, ys, ms = 1, aspect_ratio = true, legend = :none)

    plot!(pl, xs_edges, ys_edges, color = colors)

    plot!(pl, [1, 0.5, 0, 1], [0, h, 0, 0])

    return pl
end
