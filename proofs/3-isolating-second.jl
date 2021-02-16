### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 8453cee0-6dd4-11eb-247a-a7ab8e8d7592
using ArbTools, JLD, MethodOfParticularSolutions, Nemo, PaynePolygon, Plots, StaticArrays

# ╔═╡ 23d68114-6c76-11eb-373e-7b8358f16372
md"# Isolating the second eigenvalue

In this notebook we compute enclosures for the first to fourth eigenvalue, this then allows us to isolate the second eigenvalue.
"

# ╔═╡ b99db082-6dd4-11eb-2819-67e94af4dfe2
md"Load the precomputed approximations. We load both a rigorous version with using only `arb` and a non-rigorous version with `Float64` used for the plotting."

# ╔═╡ 4ef6f0ea-6dd4-11eb-094d-ebd29329adc8
domains, us, λs = let
    domains_us_λs = [
        PaynePolygon.load_eigenfunction(
            "../data/approximate-eigenfunction-$i.jld",
            T = arb,
        ) for i = 1:4
    ]

    tuple(zip(domains_us_λs...)...)
end

# ╔═╡ 27994e42-6f89-11eb-0666-f3ea843f0918
us_float64 = let
    domains_us_λs = [
        PaynePolygon.load_eigenfunction(
            "../data/approximate-eigenfunction-$i.jld",
            T = Float64,
        ) for i = 1:4
    ]

    getindex.(domains_us_λs, 2)
end

# ╔═╡ 4dfa0770-6ec8-11eb-2b22-e30892995a33
md"## Compute enclosures of eigenvalues

We need to compute enclosures of the eigenvalues for all four eigenfunctions. For this we need to bound them on the boundary and lower bound their norms. We first lower bound the norms
"

# ╔═╡ 10f2be14-6ed3-11eb-3203-b9bec4174636
md"### Lower bound the norms

We will lower bound the norm by lower bounding it on parts of the domain. For the different eigenfunctions we use different parts of the domain and in all cases we can make use of symmetries to reduce the computational cost. 

The parts are defined by piecewise linear functions. For each eigenfunction we give start and endpoints for each piece together with the area for the parts. **TODO:** Check that each area is small enough for the Faber-Krahn inequality to apply.
"

# ╔═╡ 23cbf488-6f71-11eb-071b-81096e7a1db9
A1, points1 = let i = 1
    # Parameters defining the parts of the domain
    dist = domains[i].parent(0.2)

    # Start and end point for each boundary we have to bound u on
    points = [(
        SVector(dist, zero(dist)),
        SVector(dist, tanpi(domains[i].parent(1 // 6)) * dist),
    )]

    # Area of each part
    A = 6 * dist * points[1][2][2]

    A, points
end

# ╔═╡ e09c759c-6f73-11eb-3e7f-61b9e84558fa
A2, points2 = let i = 2
    # Parameters
    dist1 = domains[i].parent(0.55)
    dist2 = domains[i].parent(0.7)
    θ₁ = domains[i].parent(1 // 6)
    θ₂ = domains[i].parent(1 // 10)
    dist3 = dist2 / cospi(θ₁ - θ₂)

    points = [
        (dist1 .* SVector(cospi(θ₁), sinpi(θ₁)), dist3 .* SVector(cospi(θ₂), sinpi(θ₂))),
        (dist3 .* SVector(cospi(θ₂), sinpi(θ₂)), dist2 .* SVector(cospi(θ₁), sinpi(θ₁))),
    ]

    A = (dist2 - dist1) * norm(points[2][2] - points[2][1])

    A, points
end

# ╔═╡ dfd06906-6f88-11eb-199b-81e3d0e9ae61
A3, points3 = let i = 3
    # Parameters
    dist1 = domains[i].parent(0.55)
    dist2 = domains[i].parent(0.7)
    θ = domains[i].parent(1 // 6)
    θ_diff = domains[i].parent(1 // 15)
    dist3 = dist2 / cospi(θ_diff)

    points = [
        (
            dist1 .* SVector(cospi(θ), sinpi(θ)),
            dist3 .* SVector(cospi(θ - θ_diff), sinpi(θ - θ_diff)),
        ),
        (
            dist3 .* SVector(cospi(θ - θ_diff), sinpi(θ - θ_diff)),
            dist3 .* SVector(cospi(θ + θ_diff), sinpi(θ + θ_diff)),
        ),
        (
            dist3 .* SVector(cospi(θ + θ_diff), sinpi(θ + θ_diff)),
            dist1 .* SVector(cospi(θ), sinpi(θ)),
        ),
    ]

    A = (dist2 - dist1) * norm(points[2][2] - points[2][1]) / 2

    A, points
end

# ╔═╡ 4ecb76b8-6f8c-11eb-0de4-3f6ca7d195b8
A4, points4 = let i = 4
    A = domains[i].parent(1) # TODO: Fix this

    dist1 = domains[i].parent(0.55)
    dist2 = domains[i].parent(0.7)
    θ = domains[i].parent(1 // 2 - 1 // 15)
    dist3 = dist2 / cospi(θ - 1 // 2)
    points = [
        (dist1 .* SVector(zero(dist1), one(dist1)), dist3 .* SVector(cospi(θ), sinpi(θ))),
        (dist3 .* SVector(cospi(θ), sinpi(θ)), dist2 .* SVector(zero(dist1), one(dist1))),
    ]

    A = (dist2 - dist1) * norm(points[2][2] - points[2][1])

    A, points
end

# ╔═╡ 30ab66b8-6f8f-11eb-3ad7-8b53de2a2cf0
begin
    As = [A1, A2, A3, A4]
    points = [points1, points2, points3, points4]
end

# ╔═╡ bb3e8878-6f73-11eb-1f6e-bd149472b122
md"In the figures below you see a (crude) approximation of the eigenfunctions with the parts where the norms are lower bounded highlighted, for symmetry reasons we only have to check the parts in red."

# ╔═╡ f2b2c20a-6ed2-11eb-2a48-15fe49e394a2
let i = 1
    pl = eigenfunctionheatmap(domains[i], us_float64[i], λs[i])

    M = j -> [cospi(j / 3) sinpi(j / 3); -sinpi(j / 3) cospi(j / 3)]
    pts = [M(j) * points1[1][2] for j = 0:6]
    plot!(
        pl,
        Float64.(getindex.(pts, 1)),
        Float64.(getindex.(pts, 2)),
        label = "",
        color = :black,
        linewidth = 2,
    )

    for (start, stop) in points1
        plot!(
            pl,
            Float64[start[1], stop[1]],
            Float64[start[2], stop[2]],
            label = "",
            color = :red,
            linewidth = 2,
        )
    end

    pl
end

# ╔═╡ eeba9d98-6f73-11eb-207a-9f91cf23bc4a
let i = 2
    pl = eigenfunctionheatmap(domains[i], us_float64[i], λs[i])

    points = Float64[points2[1][1] points2[1][2] (
        points2[1][2] + 2(points2[2][2] - points2[2][1])
    ) points2[1][1]]

    M = j -> [cospi(j / 3) sinpi(j / 3); -sinpi(j / 3) cospi(j / 3)]
    for j = 0:5
        pts = M(j) * points
        plot!(pl, pts[1, :], pts[2, :], label = "", color = :black, linewidth = 2)
    end

    for (start, stop) in points2
        plot!(
            pl,
            Float64[start[1], stop[1]],
            Float64[start[2], stop[2]],
            label = "",
            color = :red,
            linewidth = 2,
        )
    end

    pl
end

# ╔═╡ f866e5c6-6f88-11eb-1bad-a969825ae54e
let i = 3
    pl = eigenfunctionheatmap(domains[i], us_float64[i], λs[i])

    pts = [points3[1][1] points3[2][1] points3[3][1] points3[1][1]]
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]
        plot!(
            pl,
            a * Float64.(pts[1, :]),
            b * Float64.(pts[2, :]),
            label = "",
            color = :black,
            linewidth = 2,
        )
    end

    for (start, stop) in points3
        plot!(
            pl,
            Float64[start[1], stop[1]],
            Float64[start[2], stop[2]],
            label = "",
            color = :red,
            linewidth = 2,
        )
    end

    pl
end

# ╔═╡ 2e88aa42-6f8c-11eb-1b6b-f592a9760507
let i = 4
    pl = eigenfunctionheatmap(domains[i], us_float64[i], λs[i])

    pts = Float64[points4[1][1] points4[1][2] (
        points4[1][2] + 2(points4[2][2] - points4[2][1])
    ) points4[1][1]]
    for b in [1, -1]
        plot!(
            pl,
            Float64.(pts[1, :]),
            b * Float64.(pts[2, :]),
            label = "",
            color = :black,
            linewidth = 2,
        )
    end

    for (start, stop) in points4
        plot!(
            pl,
            Float64[start[1], stop[1]],
            Float64[start[2], stop[2]],
            label = "",
            color = :red,
            linewidth = 2,
        )
    end

    pl
end

# ╔═╡ 0ebf9b3c-6f8f-11eb-355b-4d7b4c851c40
md"Now we compute lower bounds of the absolute value on the parts in red"

# ╔═╡ 9ad82c38-6f8f-11eb-3441-056e38ab7675
norms2_lower = let
    norms2 = Vector{arb}(undef, length(domains))

    # Compute lower bound of u^2 on the relevant parts for each eigenfunction
    for i = 1:4
        res = domains[i].parent(Inf)
        for pts in points[i]
            v = pts[2] - pts[1]
            p(t) = pts[1] + t .* v

            res_tmp = domains[i].parent(Inf)
            N = 0
            evals = Vector{arb}(undef, N + 1)
            Threads.@threads for j = 0:N
                evals[j+1] = us[i](p(j / N), λs[i])^2
            end
            for e in evals
                res_tmp = min(res_tmp, e)
            end

            if true
                res_tmp =
                    -enclosemaximum(
                        t -> -(us[i](p(t), λs[i]))^2,
                        domains[i].parent(0),
                        domains[i].parent(1),
                        evaltype = :taylor,
                        n = length(coefficients(us[i])) ÷ 2,
                        rtol = 1e-1,
                        maxevals = 20,
                        show_trace = true,
                    )
            end
            @show i res_tmp
            res = min(res, res_tmp)
        end

        norms2[i] = res * As[i]
    end

    # Multiply with number of copies of the domain
    norms2 .*= [1, 6, 4, 2]
    #norms[2:4] .= domains[1].parent(0)

    norms2
end

# ╔═╡ efcd79ce-6f9b-11eb-1caa-933b5a5eaa33
norms = sqrt.(norms2_lower)

# ╔═╡ 3bf6a08e-6dd5-11eb-056a-79ed2b8f6bc6
norms_approx = [
    MethodOfParticularSolutions.norm(domains[i], us[i], λs[i], numpoints = 100) for
    i in eachindex(domains)
]

# ╔═╡ 08fd1908-6f95-11eb-3481-07fe6869b396
Float64.(norms_approx ./ norms)

# ╔═╡ 6147348e-6dd5-11eb-3cc3-6595d4e1fe57
md"### Upper bound the value on the boundary

For the boundary we don't need the tightest bound possible. Instead we compute an approximate bound by evaluating on a number of points and then we prove that they are bounded by a small multiple of this."

# ╔═╡ 8ca4b5c0-6dd5-11eb-105c-c3a165f91b93
approx_max = let
    compute_approx_max(domain, u, λ) = begin
        max_numpoints = length(coefficients(u)) # TODO: Increase this
        pts, bds = boundary_points(domain, u, length(coefficients(u)), max_numpoints)
        values = similar(pts, arb)
        Threads.@threads for i in eachindex(pts)
            values[i] = u(pts[i], λ, boundary = bds[i])
        end
        m = zero(λ)
        for v in values
            m = max(m, abs(v))
        end
        return m
    end
    [compute_approx_max(domains[i], us[i], λs[i]) for i in eachindex(domains)]
end

# ╔═╡ 5dffe32a-6e0c-11eb-0ac2-ff41ca816c6b
md"Not we prove that they are bounded by twice this approximate value on the boundary. For symmetry reasons we only have to bound them on a small subset of the boundary. **This is not quite what we do, we try to improve the bound until we get one which is smaller than twice the approximate value**."

# ╔═╡ bafc2126-6ec8-11eb-0ab4-23f0e59beb12
active_boundaries = [
    MethodOfParticularSolutions.active_boundaries(domains[i], us[i]) for
    i in eachindex(domains)
]

# ╔═╡ 46c2ca40-6e34-11eb-3560-5f084a97f76f
md"Below are the part of the boundary that we havea to check for the first/second and the third/fourth respectively."

# ╔═╡ 6768fefc-6e33-11eb-0e3e-6f7ead6a152c
let i = 1
    pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)
    for boundary in active_boundaries[i]
        p(t) = boundary_parameterization(t, domains[i], boundary)
        start = Float64.(p(0))
        stop = Float64.(p(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1)))
        plot!(
            pl,
            [start[1], stop[1]],
            [start[2], stop[2]],
            color = :red,
            linewidth = 3,
            title = "Boundary to check for first and second eigenvalue",
        )
    end
    pl
end

# ╔═╡ 26a68726-6e34-11eb-1d44-e5bf2700d36f
let i = 3
    pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)
    for boundary in active_boundaries[i]
        p(t) = boundary_parameterization(t, domains[i], boundary)
        start = Float64.(p(0))
        stop = Float64.(p(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1)))
        plot!(
            pl,
            [start[1], stop[1]],
            [start[2], stop[2]],
            color = :red,
            linewidth = 3,
            title = "Boundary to check for third and fourth eigenvalue",
        )
    end
    pl
end

# ╔═╡ b315e674-6e0d-11eb-0c4d-49f4bac913d9
ok₁, res₁ = let i = 1
    ok = true
    res = domains[i].parent(0)
    for boundary in active_boundaries[i]
        stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1))
        p(t) = boundary_parameterization(t, domains[i], boundary)
        ok_tmp, res_tmp = PaynePolygon.bounded_by(
            t -> us[i](p(t), λs[i]; boundary),
            domains[i].parent(0), # TODO: Change this
            domains[i].parent(stop),
            2approx_max[i],
            use_taylor = :true,
            n = length(coefficients(us[i])),
            start_intervals = 8,
            show_trace = true,
            show_evaluations = true,
            return_enclosure = true,
        )
        ok = ok && ok_tmp
        res = max(res, res_tmp)
    end
    ok, res
end

# ╔═╡ 8228d29a-6e35-11eb-3a2a-f3fa3c31e8a1
ok₂, res₂ = let i = 2
    ok = true
    res = domains[i].parent(0)
    for boundary in active_boundaries[i]
        stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1))
        p(t) = boundary_parameterization(t, domains[i], boundary)
        ok_tmp, res_tmp = PaynePolygon.bounded_by(
            t -> us[i](p(t), λs[i]; boundary),
            domains[i].parent(stop), # TODO: Change this
            domains[i].parent(stop),
            2approx_max[i],
            use_taylor = :true,
            n = length(coefficients(us[i])),
            start_intervals = 8,
            show_trace = true,
            show_evaluations = true,
            return_enclosure = true,
        )
        ok = ok && ok_tmp
        res = max(res, res_tmp)
    end
    ok, res
end

# ╔═╡ b48e4124-6e35-11eb-2321-bd2ab359d7c5
ok₃, res₃ = let i = 3
    ok = true
    res = domains[i].parent(0)
    for boundary in active_boundaries[i]
        stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1))
        p(t) = boundary_parameterization(t, domains[i], boundary)
        ok_tmp, res_tmp = PaynePolygon.bounded_by(
            t -> us[i](p(t), λs[i]; boundary),
            domains[i].parent(stop),  # TODO: Change this
            domains[i].parent(stop),
            2approx_max[i],
            use_taylor = :true,
            n = length(coefficients(us[i])),
            start_intervals = 8,
            show_trace = true,
            show_evaluations = true,
            return_enclosure = true,
        )
        ok = ok && ok_tmp
        res = max(res, res_tmp)
    end
    ok, res
end

# ╔═╡ c5d5a696-6e35-11eb-101c-4555a197b005
ok₄, res₄ = let i = 4
    ok = true
    res = domains[i].parent(0)
    for boundary in active_boundaries[i]
        stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1 // 2, 1))
        p(t) = boundary_parameterization(t, domains[i], boundary)
        ok_tmp, res_tmp = PaynePolygon.bounded_by(
            t -> us[i](p(t), λs[i]; boundary),
            domains[i].parent(stop),  # TODO: Change this
            domains[i].parent(stop),
            2approx_max[i],
            use_taylor = :true,
            n = length(coefficients(us[i])),
            start_intervals = 8,
            show_trace = true,
            show_evaluations = true,
            return_enclosure = true,
        )
        ok = ok && ok_tmp
        res = max(res, res_tmp)
    end
    ok, res
end

# ╔═╡ 15cb030c-6e37-11eb-2d26-c7676a668ed3
md"Check that all of them in fact succeded."

# ╔═╡ 0bd4cd1a-6e37-11eb-13a2-b11a437fb728
if all([ok₁, ok₂, ok₃, ok₄])
    md"All succeded!"
else
    md"Not all of them succeded!"
end

# ╔═╡ 0193502e-6e37-11eb-0959-0309ab0e8a81
max_values = [res₁, res₂, res₃, res₄]

# ╔═╡ 61ce8ef4-6e37-11eb-01fc-1f2c3988c964
md"Finally we can compute the enclosures!"

# ╔═╡ d2c27968-6f9a-11eb-1af6-8d408aa330a6
μs = [
    sqrt(MethodOfParticularSolutions.area(domains[i])) * max_values[i] / norms[i] for
    i in eachindex(domains)
]

# ╔═╡ ed377550-6f9a-11eb-365f-855367472a3f
Float64.(μs)

# ╔═╡ ffd288e0-6ec7-11eb-2fa0-058f7b36c70b
enclosures = [
    begin
        lower = λs[i] / (1 + getinterval(μs[i])[2])
        upper = λs[i] / (1 - getinterval(μs[i])[2])
        setinterval(lower, upper)
    end for i in eachindex(domains)
]

# ╔═╡ 301baa02-6f9c-11eb-00f7-af445359b0a0
Float64.(radius.(enclosures))

# ╔═╡ f598f6a6-6ec8-11eb-24d9-1d5ffd8bb2bd
md"## Handle the 2-cluster

The third and fourth eigenvalue overlaps so instead we have to follow the procedure in the lemma in the paper.
"

# ╔═╡ d4a649c0-6e38-11eb-1772-657a8c0e5b00
overlaps(enclosures[3], enclosures[4])

# ╔═╡ eb4f83b2-6e38-11eb-1654-9921d470af28
md"Let"

# ╔═╡ 981590c8-6e39-11eb-0a0a-bb76d3df6a25
Λ = setunion(enclosures[3], enclosures[4])

# ╔═╡ a661e938-6e39-11eb-2f73-bf70c7520d22
md"We will prove that there has to be at least two eigenvalues in the interval"

# ╔═╡ 02269288-6e39-11eb-0004-bfe8a0065198
Λ´ = ball(midpoint(Λ), 2radius(Λ))

# ╔═╡ 6693a8fa-6e39-11eb-07a8-c9868e805742
md"If there are not two eigenvalues in this interval then we get as a lower bound for $\alpha$"

# ╔═╡ bace2814-6e39-11eb-1d20-ad0f78c2721f
α = radius(Λ)

# ╔═╡ d43e9932-6e39-11eb-143c-65bbf3f9f84d
md"For $g(x)$ we have the bound"

# ╔═╡ d8672e3e-6e39-11eb-392f-8b6ec139f7da
g = sqrt(2MethodOfParticularSolutions.area(domains[1])) / 4domains[1].parent(π)

# ╔═╡ f1b47914-6e39-11eb-23ee-bd7dc6ef3387
md"Giving us the $L^\infty$ bounds"

# ╔═╡ 2c9d1e64-6e3a-11eb-2bf6-71e8222ee966
bound₃ = max_values[3] * (1 + g * λs[3] * (1 / (1 - μs[3]) + 1 / α * (1 + μs[3]^2 / α^2)))

# ╔═╡ 48ee4a0c-6e3a-11eb-1d20-d9939ae5ab88
bound₄ = max_values[4] * (1 + g * λs[4] * (1 / (1 - μs[4]) + 1 / α * (1 + μs[4]^2 / α^2)))

# ╔═╡ 8b978df0-6e3a-11eb-075e-072925439fe5
md"We now evaluate both `u[3]` and `u[4]` at the points"

# ╔═╡ 951e7df4-6e3a-11eb-1c26-798f0b393333
two_points = [domains[3].parent.([0.5, 0.5]), domains[3].parent.([-0.5, 0.5])]

# ╔═╡ 530d4704-6e3c-11eb-1103-c3cbe71e81d4
md"and check that we can definitely determine the signs at them"

# ╔═╡ 668266e6-6e3c-11eb-2b65-699b9eb978f6
could_determine_sign =
    abs(us[3](two_points[1], λs[3])) - bound₃ > 0 &&
    abs(us[3](two_points[2], λs[3])) - bound₃ > 0 &&
    abs(us[4](two_points[1], λs[3])) - bound₄ > 0 &&
    abs(us[4](two_points[2], λs[3])) - bound₄ > 0

# ╔═╡ e514550a-6e3c-11eb-1724-439eb0841b50
md"Finally check that $u_3$ has different signs at the two points and that $u_4$ has the same sign"

# ╔═╡ f6c3fb8e-6e3c-11eb-1df7-95194b512b20
correct_signs =
    us[3](two_points[1], λs[3]) * us[3](two_points[2], λs[3]) < 0 &&
    us[4](two_points[1], λs[4]) * us[4](two_points[2], λs[4]) > 0

# ╔═╡ dd03a248-6e3d-11eb-20b3-71fdf50a05b0
md"## Conclude

To finish we check that everything actually succeded and so that the first and second eigenvalue do not overlap with $\Lambda'$"

# ╔═╡ c077b344-6f9b-11eb-08c4-abb4033bc203
enclosures[1], enclosures[2], Λ´, Float64(ArbTools.lbound(Λ´ - λs[2]))

# ╔═╡ cdb55222-6e93-11eb-1cbf-4b9f17d72c2b
begin
    λ₁_isolated = !overlaps(enclosures[1], enclosures[2]) && !overlaps(enclosures[1], Λ´)
    λ₂_isolated = !overlaps(enclosures[2], enclosures[1]) && !overlaps(enclosures[2], Λ´)
    if λ₁_isolated && λ₂_isolated && could_determine_sign && correct_signs
        md"Everything was successful! 🎈🎈🎈"
    else
        md"The proof did not succeed 😦"
    end
end

# ╔═╡ e38d5874-6e98-11eb-0d2a-95d969e82b42
md"Finally we save the data that is used in the paper"

# ╔═╡ 1a4c99ee-6e99-11eb-170d-4d40faa1a61f
let dir = "../data"
    save(
        joinpath(dir, "enclosures.jld"),
        "enclosures_dump",
        ArbTools.arb_dump.(enclosures),
        "norms_dump",
        ArbTools.arb_dump.(norms),
        "max_boundary_dump",
        ArbTools.arb_dump.(max_values),
        "μs_dump",
        ArbTools.arb_dump.(μs),
    )
    save(
        joinpath(dir, "cluster.jld"),
        "Λ´_dump",
        ArbTools.arb_dump(Λ´),
        "L_infinity_bound_u₃_dump",
        ArbTools.arb_dump(bound₃),
        "L_infinity_bound_u₄_dump",
        ArbTools.arb_dump(bound₄),
    )
end

# ╔═╡ Cell order:
# ╟─23d68114-6c76-11eb-373e-7b8358f16372
# ╠═8453cee0-6dd4-11eb-247a-a7ab8e8d7592
# ╠═b99db082-6dd4-11eb-2819-67e94af4dfe2
# ╠═4ef6f0ea-6dd4-11eb-094d-ebd29329adc8
# ╠═27994e42-6f89-11eb-0666-f3ea843f0918
# ╟─4dfa0770-6ec8-11eb-2b22-e30892995a33
# ╟─10f2be14-6ed3-11eb-3203-b9bec4174636
# ╠═23cbf488-6f71-11eb-071b-81096e7a1db9
# ╟─e09c759c-6f73-11eb-3e7f-61b9e84558fa
# ╟─dfd06906-6f88-11eb-199b-81e3d0e9ae61
# ╟─4ecb76b8-6f8c-11eb-0de4-3f6ca7d195b8
# ╠═30ab66b8-6f8f-11eb-3ad7-8b53de2a2cf0
# ╟─bb3e8878-6f73-11eb-1f6e-bd149472b122
# ╟─f2b2c20a-6ed2-11eb-2a48-15fe49e394a2
# ╟─eeba9d98-6f73-11eb-207a-9f91cf23bc4a
# ╟─f866e5c6-6f88-11eb-1bad-a969825ae54e
# ╟─2e88aa42-6f8c-11eb-1b6b-f592a9760507
# ╟─0ebf9b3c-6f8f-11eb-355b-4d7b4c851c40
# ╠═9ad82c38-6f8f-11eb-3441-056e38ab7675
# ╠═efcd79ce-6f9b-11eb-1caa-933b5a5eaa33
# ╠═3bf6a08e-6dd5-11eb-056a-79ed2b8f6bc6
# ╠═08fd1908-6f95-11eb-3481-07fe6869b396
# ╟─6147348e-6dd5-11eb-3cc3-6595d4e1fe57
# ╠═8ca4b5c0-6dd5-11eb-105c-c3a165f91b93
# ╟─5dffe32a-6e0c-11eb-0ac2-ff41ca816c6b
# ╠═bafc2126-6ec8-11eb-0ab4-23f0e59beb12
# ╟─46c2ca40-6e34-11eb-3560-5f084a97f76f
# ╟─6768fefc-6e33-11eb-0e3e-6f7ead6a152c
# ╟─26a68726-6e34-11eb-1d44-e5bf2700d36f
# ╠═b315e674-6e0d-11eb-0c4d-49f4bac913d9
# ╠═8228d29a-6e35-11eb-3a2a-f3fa3c31e8a1
# ╠═b48e4124-6e35-11eb-2321-bd2ab359d7c5
# ╠═c5d5a696-6e35-11eb-101c-4555a197b005
# ╟─15cb030c-6e37-11eb-2d26-c7676a668ed3
# ╠═0bd4cd1a-6e37-11eb-13a2-b11a437fb728
# ╠═0193502e-6e37-11eb-0959-0309ab0e8a81
# ╟─61ce8ef4-6e37-11eb-01fc-1f2c3988c964
# ╠═d2c27968-6f9a-11eb-1af6-8d408aa330a6
# ╠═ed377550-6f9a-11eb-365f-855367472a3f
# ╠═301baa02-6f9c-11eb-00f7-af445359b0a0
# ╠═ffd288e0-6ec7-11eb-2fa0-058f7b36c70b
# ╟─f598f6a6-6ec8-11eb-24d9-1d5ffd8bb2bd
# ╠═d4a649c0-6e38-11eb-1772-657a8c0e5b00
# ╟─eb4f83b2-6e38-11eb-1654-9921d470af28
# ╠═981590c8-6e39-11eb-0a0a-bb76d3df6a25
# ╟─a661e938-6e39-11eb-2f73-bf70c7520d22
# ╠═02269288-6e39-11eb-0004-bfe8a0065198
# ╟─6693a8fa-6e39-11eb-07a8-c9868e805742
# ╠═bace2814-6e39-11eb-1d20-ad0f78c2721f
# ╟─d43e9932-6e39-11eb-143c-65bbf3f9f84d
# ╠═d8672e3e-6e39-11eb-392f-8b6ec139f7da
# ╟─f1b47914-6e39-11eb-23ee-bd7dc6ef3387
# ╠═2c9d1e64-6e3a-11eb-2bf6-71e8222ee966
# ╠═48ee4a0c-6e3a-11eb-1d20-d9939ae5ab88
# ╟─8b978df0-6e3a-11eb-075e-072925439fe5
# ╠═951e7df4-6e3a-11eb-1c26-798f0b393333
# ╟─530d4704-6e3c-11eb-1103-c3cbe71e81d4
# ╠═668266e6-6e3c-11eb-2b65-699b9eb978f6
# ╟─e514550a-6e3c-11eb-1724-439eb0841b50
# ╠═f6c3fb8e-6e3c-11eb-1df7-95194b512b20
# ╟─dd03a248-6e3d-11eb-20b3-71fdf50a05b0
# ╠═c077b344-6f9b-11eb-08c4-abb4033bc203
# ╠═cdb55222-6e93-11eb-1cbf-4b9f17d72c2b
# ╟─e38d5874-6e98-11eb-0d2a-95d969e82b42
# ╠═1a4c99ee-6e99-11eb-170d-4d40faa1a61f
