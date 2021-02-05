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

"""
    plot_mesh(N = 9, d, h; plot_boundary = true, highlight_fundamental = plot_boundary)

Plot the mesh given by the triple `(N, d, h)`, the same that is used
in [`stiffness_matrix`](@ref).

If `plot_boundary` is true then plot the boundary of the domain in
addition to the edges in the mesh. If `highlight_fundamental` is true
then highlight the fundamental domain to make it easier to visualize.
If `zoom_hole` is true then zoom in on one of the holes to make it
easier to see the details. If `plot_mesh` is false then don't plot the
mesh, this is useful for plotting only the boundary of the domain.
"""
function plot_mesh(
    N::Integer = 9,
    d::Integer = 0,
    h::Integer = 0;
    plot_mesh = true,
    plot_boundary = true,
    highlight_fundamental = false,
    zoom_hole = false,
)
    num_triangles = N^2
    xs, ys = mesh_triangle_midpoints(N)
    edges = collect(zip(edges_to_triangles(N)...))
    boundary, interior, vertices = boundary_and_interior(N, d, h)

    filter!(
        ((t1, t2),) ->
            !((t1 % num_triangles ∈ interior) || (t2 % num_triangles ∈ interior)) &&
                !((t1 % num_triangles ∈ vertices) && (t2 % num_triangles ∈ vertices)),
        edges,
    )

    xs_edges = zeros(2, length(edges))
    ys_edges = zeros(2, length(edges))
    colors = Array{Symbol}(undef, 1, length(edges))

    for (i, (t1, t2)) in enumerate(edges)
        if (t1 % num_triangles ∈ boundary) && (t2 % num_triangles ∈ boundary)
            # The edge touches two triangles on the boundary. We want
            # to print only half of the edge, but we have to determine
            # which half.
            e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]])

            # Determine which part of the edge we want to keep but
            # estimating which row it's on with sqrt(t1 %
            # num_triangles).
            if sqrt(t1 % num_triangles) < d + h ÷ 3
                e2 = (e1 + e2) / 2
            else
                e1 = (e1 + e2) / 2
            end
            colors[i] = :red
        elseif (t1 % num_triangles ∈ boundary) || (t2 % num_triangles ∈ boundary)
            # The edge touches one triangle on the boundary
            e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]])
            colors[i] = :red
        else
            e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]])
            colors[i] = :blue
        end
        xs_edges[1, i] = e1[1]
        xs_edges[2, i] = e2[1]
        ys_edges[1, i] = e1[2]
        ys_edges[2, i] = e2[2]
    end

    if plot_mesh
        pl = plot(xs_edges, ys_edges, color = colors, aspect_ratio = true, legend = :none)
    else
        pl = plot(aspect_ratio = true, legend = :none)
    end

    if zoom_hole
        xlims!(pl, ((d - 2) / N, (d + h + 2) / N))
        ylims!(pl, sqrt(3) / 2 .* (-(h + 1) / N, (h + 1) / N))
    end

    H = sqrt(3) / 2
    if plot_boundary
        # Boundary of hexagon
        plot!(pl, [1, 0.5, -0.5, -1, -0.5, 0.5, 1], [0, H, H, 0, -H, -H, 0], c = :black)

        # Boundary of interior triangles
        pts = [d/N -H*2h/3N; (d+h)/N 0; d/N H*2h/3N; d/N -H*2h/3N]'
        for i = 0:5
            # Rotate points by i*π/3
            M = [cospi(i / 3) sinpi(i / 3); -sinpi(i / 3) cospi(i / 3)]
            pts_rotated = M * pts
            plot!(pl, pts_rotated[1, :], pts_rotated[2, :], c = :black)
        end
    end
    if highlight_fundamental
        plot!(pl, [1, 0.5, 0, 1], [0, H, 0, 0], c = :green)
    end

    return pl
end
