### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 5b800cd2-6c55-11eb-0ccc-25420640362e
using Revise

# ╔═╡ 5dbb9e30-6c55-11eb-36a3-210e1d52a1b5
using ArbTools, JLD, MethodOfParticularSolutions, Nemo, PaynePolygon, Plots, StaticArrays

# ╔═╡ 00ab445c-6c55-11eb-380f-2b64b223938e
md"# Isolating nodal line - proof

This notebook contains the code to prove that the nodal line for $u_2$ is isolated. It uses parameters computed in `isolate-nodal-line-exploration.jl`.
"

# ╔═╡ 8039961a-6c55-11eb-2072-536d9ddf52db
md"We load the precomputed approximate eigenfunction"

# ╔═╡ 8784e15e-6c55-11eb-234c-eb1f0f423bae
domain, u, λ =
    PaynePolygon.load_eigenfunction("../data/approximate-eigenfunction-2.jld", T = arb)

# ╔═╡ 889a9cd0-6c55-11eb-2a2f-09e6225e11d4
md"At the moment the code assumes the eigenfunction is positive at the center"

# ╔═╡ 96c1da6e-6c55-11eb-3969-090f05e7647a
@assert u([domain.parent(1e-8), domain.parent(0)], λ) > 0

# ╔═╡ 98472984-6c55-11eb-2f45-83a0977ca89a
distance = domain.parent(load("../data/distance.jld")["distance"])

# ╔═╡ 3b87f128-6c56-11eb-304c-7b1266e73c5c
md"Define a parameterization of the line where we want to prove that `u` is negative"

# ╔═╡ 57227a8e-6c56-11eb-3681-eb0a5b4e775d
p = let
    y_max = distance * tanpi(domain.parent(1 // 6))
    p(t::arb) = SVector{2,arb}(distance, t * y_max)
	# Hack to get the output of the correct type
	p(t::arb_series) = SVector{2,arb_series}(distance + 0t, t * y_max)
end

# ╔═╡ 2cb946a8-6c57-11eb-10e5-73e21c21b614
md"In practice the lower bound is attained at $p(0)$, use that to speed up computations later"

# ╔═╡ 39278e10-6c57-11eb-2fa5-ab37f8c00a2a
lower_bound = u(p(domain.parent(0)), λ)

# ╔═╡ 51a735ee-6c57-11eb-1ff0-17b875452545
md"Compute an enclosure of the maximum value of `u(p(t), λ)` for `0 ≤ t ≤ 1`. Notice that `u` is negative on the line so this will tell us how far away from zero it is."

# ╔═╡ a9ce8d12-6c57-11eb-0260-9d59be9e0d95
value = enclosemaximum(t -> u(p(t), λ), domain.parent(0), domain.parent(1), evaltype = :taylor, n = length(coefficients(u))÷2, rtol = 1e-2, show_trace = true; lower_bound)

# ╔═╡ 3427b152-6c58-11eb-073d-3b337dbf73a7
md"Once we have this value we can determine what $L^\infty$ bound we have to beat. Instead of computing the a as tight bound as possible we can then compute until we know the bound is good enough."

# ╔═╡ 7e4e1bd6-6c58-11eb-2ad1-958d24807648
# TODO: Compute what bound we have to beat

# ╔═╡ 84a91dd0-6c58-11eb-224f-09806f358f5f
# TODO: Prove that we satsify this bound

# ╔═╡ Cell order:
# ╟─00ab445c-6c55-11eb-380f-2b64b223938e
# ╠═5b800cd2-6c55-11eb-0ccc-25420640362e
# ╠═5dbb9e30-6c55-11eb-36a3-210e1d52a1b5
# ╟─8039961a-6c55-11eb-2072-536d9ddf52db
# ╠═8784e15e-6c55-11eb-234c-eb1f0f423bae
# ╟─889a9cd0-6c55-11eb-2a2f-09e6225e11d4
# ╠═96c1da6e-6c55-11eb-3969-090f05e7647a
# ╠═98472984-6c55-11eb-2f45-83a0977ca89a
# ╟─3b87f128-6c56-11eb-304c-7b1266e73c5c
# ╠═57227a8e-6c56-11eb-3681-eb0a5b4e775d
# ╟─2cb946a8-6c57-11eb-10e5-73e21c21b614
# ╠═39278e10-6c57-11eb-2fa5-ab37f8c00a2a
# ╟─51a735ee-6c57-11eb-1ff0-17b875452545
# ╠═a9ce8d12-6c57-11eb-0260-9d59be9e0d95
# ╟─3427b152-6c58-11eb-073d-3b337dbf73a7
# ╠═7e4e1bd6-6c58-11eb-2ad1-958d24807648
# ╠═84a91dd0-6c58-11eb-224f-09806f358f5f
