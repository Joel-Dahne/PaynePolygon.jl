### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 5791440e-6bb7-11eb-25da-31e6ec046e15
using Revise

# ╔═╡ 3d88a638-6bb7-11eb-097c-c38d8e633876
using ArbTools, JLD,
    MethodOfParticularSolutions, Nemo, NLsolve, Optim, PaynePolygon, Plots, StaticArrays

# ╔═╡ 88875ad6-6bb1-11eb-3a51-9535b2be326c
md"# Isolating the nodal line - exploration

This notebook contains the first steps to proving that the nodal line is isolated. Non of the computations here are rigorous, instead this notebook is used to setting up the problem and determining good values for the parameters in the proof. The actuall proof is carried out in a different notebook, but uses parameters computed here.
"

# ╔═╡ 2bc917de-6bb7-11eb-2c68-037d72cc58da
md"We load the precomputed approximate eigenfunction"

# ╔═╡ 37b67078-6bb7-11eb-1254-ddb708863048
domain, u, λ =
    PaynePolygon.load_eigenfunction("../data/approximate-eigenfunction-2.jld", T = Float64)

# ╔═╡ b325c45a-6c40-11eb-2408-a54a3de4a998
md"At the moment the code assumes that the eigenfunction is positive at the center."

# ╔═╡ be92b528-6c40-11eb-19d0-79a24606712b
@assert u([domain.parent(1e-8), domain.parent(0)], λ) > 0

# ╔═╡ 6eed4472-6bb7-11eb-255e-fb33f7c0f05d
md"To visualise the nodal line we plot the absolute value of the eigenfunction, zoomed in on the center of the domain. **TODO:** Make a better plot of this."

# ╔═╡ 76d8e484-6bb7-11eb-0363-457a755ddadc
pl = let num_points = 50
    eigenfunctionheatmap(
        domain,
        u,
        λ,
        range(-0.5, 0.5, length = num_points),
        range(-0.5, 0.5, length = num_points),
        absolute_value = true,
    )
end

# ╔═╡ ad5bbcc0-6bb7-11eb-2168-2b7d5e2ea48c
md"The first step to prove that the nodal line is isolated is to prove that the eigenfunction is bounded away from zero on a loop around the center. We define the loop by a straight line with a given distance from the origin. The distance is defined furher down."

# ╔═╡ 7c16c964-6bb9-11eb-3d85-55b7a29c6c91
md"The line is then extended by symmetry"

# ╔═╡ 84b13d60-6bba-11eb-01f8-318885e2a6d0
md"Due to the symmetries of the eigenfunction it's enough to prove that it's bounded away from zero on the red part of the loop. Zooming in on that part we get"

# ╔═╡ 48b55570-6bbb-11eb-3ce9-6702a18d56ed
md"We define a parameterization, $p(t)$, of the line"

# ╔═╡ 86198916-6bbb-11eb-05e7-8b45e98f08a9
md"We can the plot an approximation of the `u` along the line"

# ╔═╡ 52d22d12-6bca-11eb-1a10-adc4c2523c39
md"Zooming in closer to `t = 0` we see that it does look like it's positive."

# ╔═╡ 029aa43a-6c3f-11eb-1858-fb78b6eff99d
md"We want the value of `distance` to be such that we are as far away from zero as possible in the above plot. Since the it's the closest to zero at `(distance, 0)` we can plot the value of `u` along the line `(x, 0)`"

# ╔═╡ 4d6ae718-6c3f-11eb-3dad-e5a50eef4e86
let
    distances = range(0, 11 / 27, length = 1000)
    s(distance) = SVector(domain.parent(distance), domain.parent(0))
    plot(distances, d -> Float64(u(s(d), λ)), legend = :none)
end

# ╔═╡ db3999ae-6c3f-11eb-0bbe-e1ff1d1862ba
md"We want to find the **minimum** value for this function. In this case we don't need a rigorous value, any good approximation is fine."

# ╔═╡ 28abf718-6c40-11eb-2d2c-afe5ce807dd3
distance, value = let
    f(distance) = Float64(u(SVector(domain.parent(distance), domain.parent(0)), λ))
    res = optimize(f, 0, 11 / 27)
    res.minimizer, res.minimum
end

# ╔═╡ 707e4dd0-6bb8-11eb-3998-35d9fc524bc5
distance

# ╔═╡ 37450842-6bb8-11eb-11c2-f51cd780a816
start, stop = SVector(distance, 0), SVector(distance, tan(π / 6) * distance)

# ╔═╡ bb7a793a-6bb8-11eb-25e3-3318fbac0edf
let num_points = 50
    pl = eigenfunctionheatmap(
        domain,
        u,
        λ,
        range(-0.5, 0.5, length = num_points),
        range(-0.5, 0.5, length = num_points),
        absolute_value = true,
    )

    plot!(
        pl,
        [start[1], stop[1]],
        [start[2], stop[2]],
        label = "",
        color = :red,
        linewidth = 2,
    )
end

# ╔═╡ 85f2257a-6bb9-11eb-2432-010b6e7ebed5
let num_points = 50
    pl = eigenfunctionheatmap(
        domain,
        u,
        λ,
        range(-0.5, 0.5, length = num_points),
        range(-0.5, 0.5, length = num_points),
        absolute_value = true,
    )

    M = i -> [cospi(i / 3) sinpi(i / 3); -sinpi(i / 3) cospi(i / 3)]
    pts = [M(i) * stop for i = 0:6]
    plot!(
        pl,
        getindex.(pts, 1),
        getindex.(pts, 2),
        label = "",
        color = :black,
        linewidth = 2,
    )

    plot!(
        pl,
        [start[1], stop[1]],
        [start[2], stop[2]],
        label = "",
        color = :red,
        linewidth = 2,
    )
end

# ╔═╡ c5ebbbd4-6bba-11eb-0808-2b62067bb188
let num_points = 50
    pl = eigenfunctionheatmap(
        domain,
        u,
        λ,
        range(distance - 0.1, distance + 0.1, length = num_points),
        range(-0.1, stop[2] + 0.1, length = num_points),
    )

    plot!(
        pl,
        [start[1], stop[1]],
        [start[2], stop[2]],
        label = "",
        color = :red,
        linewidth = 2,
    )
end

# ╔═╡ 50be8ce6-6bbb-11eb-2936-69ebc659316d
p = let distance = domain.parent(distance)
    y_max = distance * tanpi(domain.parent(1 // 6))
    p(t) = SVector(distance, t * y_max)
end

# ╔═╡ a5cace5c-6bbb-11eb-00fc-75c2806143a9
let
    ts = range(0, 1, length = 1000)
    plot(ts, t -> Float64(u(p(t), λ)), legend = :none)
end

# ╔═╡ 26108832-6bca-11eb-25fd-5b82e987ea81
let
    ts = range(0, 0.1, length = 1000)
    plot(ts, t -> Float64(u(p(t), λ)), legend = :none)
    hline!([0])
end

# ╔═╡ b9c346c4-6c55-11eb-23e9-25939f697bd0
md"The value of `distance` will be used in the proof so store it for later use"

# ╔═╡ dbf1c536-6c55-11eb-3365-cd930ba42cf6
save("../data/distance.jld", "distance", distance)

# ╔═╡ effc8eb6-6c42-11eb-156b-513f9a6ea67b
md"Finally we check if it seems like the $L^\infty$ bound for `u` from Theorem 5.2 in the paper is good enough. We do this by computing an approximation of the bound and compare that with `value`, if the bound is smaller than we are good to go! The first step is to compute approximations of the the norm `n` of `u` and its maximum value on the boundary `m`."

# ╔═╡ 2917dbc4-6c43-11eb-3edc-d7d8637dca50
n, m = let
    _, n, m = MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domain,
        u,
        λ,
        store_trace = true,
        max_numpoints = 32 * length(coefficients(u)),
    )
    Float64(n), Float64(m)
end

# ╔═╡ ef47784a-6c43-11eb-0085-217b07d70666
md"With this we can compute $\mu$"

# ╔═╡ f4c21d02-6c43-11eb-10a0-7547b45f9f77
μ = Float64(sqrt(MethodOfParticularSolutions.area(domain))) * m / n

# ╔═╡ af872b10-6c43-11eb-0653-d972ceb73fc9
md"We also need a lower bound for $α$, **TODO** use a rigorously computed value, and an upper bound for $g(x)$"

# ╔═╡ c7a71caa-6c43-11eb-39b9-094b2c98efa2
α = Float64(min(λ - 31.0432, 63.7259 - λ))

# ╔═╡ 38047830-6c44-11eb-1c06-1f3382cf0960
g = sqrt(2) / π

# ╔═╡ 4c879968-6c44-11eb-322c-53c3226a9d28
md"We now get the $L^\infty$ bound as"

# ╔═╡ 56e0cb32-6c44-11eb-083c-615b69497641
bound = m * (1 + g * Float64(λ) * (1 / (1 - μ) + 1 / α * (1 + μ^2 / α^2)))

# ╔═╡ f74ffb2e-6c44-11eb-2c78-779338c92c15
md"We can now check if the bound seems to be good enough"

# ╔═╡ 11c74ed0-6c45-11eb-271e-45d7462c814a
if bound < abs(value)
    md"It seems like we are good to go!"
else
    md"The bound is to bad!"
end

# ╔═╡ 9b87666e-6c45-11eb-342f-130d021d7c20
md"If the bound is to bad we can try to estimate how small `m` would need to be for it to work."

# ╔═╡ ae70b348-6c45-11eb-3429-a18b528a9b8b
bound_from_max(m) = begin
    μ = Float64(sqrt(MethodOfParticularSolutions.area(domain))) * m / n
    return m * (1 + g * Float64(λ) * (1 / (1 - μ) + 1 / α * (1 + μ^2 / α^2)))
end

# ╔═╡ cba20082-6c45-11eb-1f9a-0bc6b246edb6
let
    pl = plot(
        10 .^ range(-10, -1, length = 100),
        bound_from_max,
        xaxis = :log10,
        yaxis = :log10,
        label = "bound_from_max(m)",
        xlabel = "m",
        legend = :topleft,
    )
    scatter!(pl, [m], [bound], label = "current bound")
    hline!(pl, [abs(value)], label = "bound to beat")
end

# ╔═╡ 99be7d3e-6c47-11eb-3ba1-335b85d1a7f3
md"The maximum that we have to get below is thus"

# ╔═╡ 50474f50-6c47-11eb-10f1-6ba469125635
m_to_beat =
    nlsolve(m -> [bound_from_max(only(m)) - abs(value)], [m], autodiff = :forward).zero

# ╔═╡ aadb0f88-6c47-11eb-20f6-f17f79b28a49
md"Or if we normalise by `n`"

# ╔═╡ b69b4270-6c47-11eb-273d-eb5885e83403
m_to_beat / n

# ╔═╡ c91220f4-6c47-11eb-2eab-01b51f39f4fe
md"Compared to the current value"

# ╔═╡ cf8cf09e-6c47-11eb-0091-870de939dd96
m / n

# ╔═╡ Cell order:
# ╟─88875ad6-6bb1-11eb-3a51-9535b2be326c
# ╠═5791440e-6bb7-11eb-25da-31e6ec046e15
# ╠═3d88a638-6bb7-11eb-097c-c38d8e633876
# ╟─2bc917de-6bb7-11eb-2c68-037d72cc58da
# ╠═37b67078-6bb7-11eb-1254-ddb708863048
# ╟─b325c45a-6c40-11eb-2408-a54a3de4a998
# ╠═be92b528-6c40-11eb-19d0-79a24606712b
# ╟─6eed4472-6bb7-11eb-255e-fb33f7c0f05d
# ╠═76d8e484-6bb7-11eb-0363-457a755ddadc
# ╟─ad5bbcc0-6bb7-11eb-2168-2b7d5e2ea48c
# ╠═707e4dd0-6bb8-11eb-3998-35d9fc524bc5
# ╠═37450842-6bb8-11eb-11c2-f51cd780a816
# ╟─bb7a793a-6bb8-11eb-25e3-3318fbac0edf
# ╟─7c16c964-6bb9-11eb-3d85-55b7a29c6c91
# ╟─85f2257a-6bb9-11eb-2432-010b6e7ebed5
# ╟─84b13d60-6bba-11eb-01f8-318885e2a6d0
# ╟─c5ebbbd4-6bba-11eb-0808-2b62067bb188
# ╟─48b55570-6bbb-11eb-3ce9-6702a18d56ed
# ╠═50be8ce6-6bbb-11eb-2936-69ebc659316d
# ╟─86198916-6bbb-11eb-05e7-8b45e98f08a9
# ╠═a5cace5c-6bbb-11eb-00fc-75c2806143a9
# ╟─52d22d12-6bca-11eb-1a10-adc4c2523c39
# ╠═26108832-6bca-11eb-25fd-5b82e987ea81
# ╟─029aa43a-6c3f-11eb-1858-fb78b6eff99d
# ╠═4d6ae718-6c3f-11eb-3dad-e5a50eef4e86
# ╟─db3999ae-6c3f-11eb-0bbe-e1ff1d1862ba
# ╠═28abf718-6c40-11eb-2d2c-afe5ce807dd3
# ╠═b9c346c4-6c55-11eb-23e9-25939f697bd0
# ╠═dbf1c536-6c55-11eb-3365-cd930ba42cf6
# ╟─effc8eb6-6c42-11eb-156b-513f9a6ea67b
# ╠═2917dbc4-6c43-11eb-3edc-d7d8637dca50
# ╟─ef47784a-6c43-11eb-0085-217b07d70666
# ╠═f4c21d02-6c43-11eb-10a0-7547b45f9f77
# ╟─af872b10-6c43-11eb-0653-d972ceb73fc9
# ╠═c7a71caa-6c43-11eb-39b9-094b2c98efa2
# ╠═38047830-6c44-11eb-1c06-1f3382cf0960
# ╟─4c879968-6c44-11eb-322c-53c3226a9d28
# ╠═56e0cb32-6c44-11eb-083c-615b69497641
# ╟─f74ffb2e-6c44-11eb-2c78-779338c92c15
# ╠═11c74ed0-6c45-11eb-271e-45d7462c814a
# ╟─9b87666e-6c45-11eb-342f-130d021d7c20
# ╠═ae70b348-6c45-11eb-3429-a18b528a9b8b
# ╠═cba20082-6c45-11eb-1f9a-0bc6b246edb6
# ╟─99be7d3e-6c47-11eb-3ba1-335b85d1a7f3
# ╠═50474f50-6c47-11eb-10f1-6ba469125635
# ╟─aadb0f88-6c47-11eb-20f6-f17f79b28a49
# ╠═b69b4270-6c47-11eb-273d-eb5885e83403
# ╟─c91220f4-6c47-11eb-2eab-01b51f39f4fe
# ╠═cf8cf09e-6c47-11eb-0091-870de939dd96
