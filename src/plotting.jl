"""
    mesh_triangle_midpoints(N::Integer = 36)

Return the midpoints of the triangles in the mesh.
"""
function mesh_triangle_midpoints(N::Integer = 36)
    num_triangles = N^2

    xs = zeros(num_triangles)
    ys = zeros(num_triangles)

    # Height and side length of each triangle
    height = sqrt(3)/2N
    side = 1/N

    # x and y increment along each row
    x_inc = cospi(2//3)*side/2
    y_inc = sinpi(2//3)*side/2

    for row in 1:N
        row_start = (row - 1)^2

        # First part of the triangles

        # x and y coordinates for the first triangle in the row
        x₀ = side*(row - 1) + side/2
        y₀ = side*tan(π/6)/2

        for i in 1:2:(2row - 1)
            xs[row_start + i] = x₀ + (i - 1)*x_inc
            ys[row_start + i] = y₀ + (i - 1)*y_inc
        end

        # Second part of the triangles

        # x and y coordinates for the second triangle in the row
        x₀ = side*(row - 1)
        y₀ = height - side*tan(π/6)/2

        for i in 2:2:(2row - 1)
            xs[row_start + i] = x₀ + (i - 2)*x_inc
            ys[row_start + i] = y₀ + (i - 2)*y_inc
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

    M = [cospi(1//3) -sinpi(1//3); sinpi(1//3) cospi(1//3)]

    # c1 + (v rotated π/3 positive direction)
    e1 = c1 + M*v
    # c1 + (v rotated -π/3 positive direction)
    e2 = c1 + inv(M)*v

    return e1, e2
end

function plot_mesh(N::Integer)
    xs, ys = mesh_triangle_midpoints(N)
    edge_to_triangle1, edge_to_triangle2 = edges_to_triangles(N)

    pl = plot([0, 1, 0.5, 0], [0, 0, sqrt(3)/2, 0], aspect_ratio = true, legend = :none)
    scatter!(pl, xs, ys)

    num_edges = div(3N^2 - 3N, 2)
    xs_edges = zeros(2, num_edges)
    ys_edges = zeros(2, num_edges)
    i = 1

    for (t1, t2) in zip(edge_to_triangle1, edge_to_triangle2)
        if t1 <= N^2 && t2 <= N^2
            e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]]);
            xs_edges[1, i] = e1[1]
            xs_edges[2, i] = e2[1]
            ys_edges[1, i] = e1[2]
            ys_edges[2, i] = e2[2]
            i += 1
        end
    end

    plot!(pl, xs_edges, ys_edges, color = :red)

    return pl
end
