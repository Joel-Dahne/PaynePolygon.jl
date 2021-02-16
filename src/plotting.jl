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

If `optimize` is true then don't plot each boundary segment separately
but plot them as longer lines. This reduces the space required for the
resulting plot substantially, in particular for finer meshes.
"""
function plot_mesh(
    N::Integer = 9,
    d::Integer = 0,
    h::Integer = 0;
    plot_mesh = true,
    plot_boundary = true,
    highlight_fundamental = false,
    zoom_hole = false,
    optimize = true,
)
    pl = plot(aspect_ratio = true, legend = :none)

    if plot_mesh
        num_triangles = N^2
        xs, ys = mesh_triangle_midpoints(N)
        edges = collect(zip(edges_to_triangles(N)...))
        boundary, interior, vertices = boundary_and_interior(N, d, h)

        if optimize
            H = sqrt(3) / 2

            ### Plot blue lines
            # Vertical blue lines
            ys_lines = let ys = range(-H, H, length = 2N + 1)[2:end-1]
                permutedims([ys ys])
            end
            xs_lines =
                [ifelse(y > 0, y / sqrt(3) - 1, -y / sqrt(3) - 1) for y in ys_lines] .* [1, -1]
            plot!(xs_lines, ys_lines, color = :blue)

            # Diagonal blue lines
            xs_lines =
                let xs = [
                        range(-1, -0.5, length = N + 1)[2:end]
                        range(-0.5, 0.5, length = N + 1)[2:end-1]
                    ]
                    permutedims([xs -reverse(xs)])
                end
            ys_lines = let ys = [range(0, H, length = N + 1)[2:end]; fill(H, N - 1)]
                permutedims([ys -reverse(ys)])
            end
            plot!(xs_lines, ys_lines, color = :blue)
            plot!(xs_lines, -ys_lines, color = :blue)

            # Red lines
            filter!(
                ((t1, t2),) ->
                    !((t1 % num_triangles ∈ interior) || (t2 % num_triangles ∈ interior)) &&
                        !(
                            (t1 % num_triangles ∈ vertices) &&
                            (t2 % num_triangles ∈ vertices)
                        ) &&
                        (
                            (t1 % num_triangles ∈ boundary) ||
                            (t2 % num_triangles ∈ boundary)
                        ),
                edges,
            )

            xs_edges = zeros(2, length(edges))
            ys_edges = zeros(2, length(edges))

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
                else
                    (t1 % num_triangles ∈ boundary) || (t2 % num_triangles ∈ boundary)
                    # The edge touches one triangle on the boundary
                    e1, e2 = PaynePolygon.triangle_edge([xs[t1], ys[t1]], [xs[t2], ys[t2]])
                end
                xs_edges[1, i] = e1[1]
                xs_edges[2, i] = e2[1]
                ys_edges[1, i] = e1[2]
                ys_edges[2, i] = e2[2]
            end

            plot!(pl, xs_edges, ys_edges, color = :red, aspect_ratio = true, legend = :none)

            # Handle interior of holes
            pts = [d/N -H*2h/3N; (d+h)/N 0; d/N H*2h/3N; d/N -H*2h/3N]'
            xs_lines = begin
                xs = range(0, h / N, length = 2h ÷ 3 + 1)[2:end-1]
                xs_lines = permutedims(
                    [fill(d / N, 4h ÷ 3 - 1) (d / N .+ [xs; h / N; reverse(xs)])],
                )

                ys = range(-H * 2h / 3N, H * 2h / 3N, length = 4h ÷ 3 + 1)[2:end-1]
                ys_lines = permutedims([ys ys])
            end
        else
            filter!(
                ((t1, t2),) ->
                    !((t1 % num_triangles ∈ interior) || (t2 % num_triangles ∈ interior)) &&
                        !(
                            (t1 % num_triangles ∈ vertices) &&
                            (t2 % num_triangles ∈ vertices)
                        ),
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

            plot!(
                pl,
                xs_edges,
                ys_edges,
                color = colors,
                aspect_ratio = true,
                legend = :none,
            )
        end
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
            if optimize
                plot!(
                    pl,
                    pts_rotated[1, :],
                    pts_rotated[2, :],
                    c = :black,
                    fill = (0, :white),
                )
            else
                plot!(pl, pts_rotated[1, :], pts_rotated[2, :], c = :black)
            end
        end
    end
    if highlight_fundamental
        plot!(pl, [1, 0.5, 0, 1], [0, H, 0, 0], c = :green)
    end

    return pl
end

plot_eigenfunction(
    domain,
    u,
    λ,
    num_xs::Integer = 50,
    num_ys::Integer = 50;
    kwargs...
) = plot_eigenfunction(
    domain,
    u, λ, range(-1, 1, length = num_xs), range(-sqrt(3) / 2, sqrt(3) / 2, length = num_ys); kwargs...)

function plot_eigenfunction(
    domain,
    u,
    λ,
    xs::AbstractVector,
    ys::AbstractVector;
    twosided = true,
    seriescolor = ifelse(twosided, :delta, :viridis)
)
    pts = SVector.(domain.parent.(xs'), domain.parent.(ys))
    res = similar(pts, Float64)
    let λ = domain.parent(λ)
        Threads.@threads for i in eachindex(pts)
            if pts[i] ∈ domain
                res[i] = u(pts[i], λ)
            else
                res[i] = NaN
            end
        end
    end

    if twosided
        m = maximum(abs, filter(!isnan, res))
        clims = (-m, m)
    else
        clims = (NaN, NaN)
    end
    pl = plot(aspect_ratio = true, legend = :none, axis = ([], false); clims)

    heatmap!(pl, xs, ys, res; seriescolor)

    H = sqrt(3) / 2
    # Boundary of hexagon
    plot!(pl, [1, 0.5, -0.5, -1, -0.5, 0.5, 1], [0, H, H, 0, -H, -H, 0], color = :black, linewidth = 2)

    # Boundary of interior triangles
    pts = let N = 27, d = 11, h = 6
        [d/N -H*2h/3N; (d+h)/N 0; d/N H*2h/3N; d/N -H*2h/3N]'
    end
    for i = 0:5
        # Rotate points by i*π/3
        M = [cospi(i / 3) sinpi(i / 3); -sinpi(i / 3) cospi(i / 3)]
        pts_rotated = M * pts
        plot!(pl, pts_rotated[1, :], pts_rotated[2, :], color = :black, linewidth = 2)
    end

    return pl
end
