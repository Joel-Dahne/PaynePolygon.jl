### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 5dbb9e30-6c55-11eb-36a3-210e1d52a1b5
using ArbTools, JLD, MethodOfParticularSolutions, Nemo, PaynePolygon, StaticArrays

# ╔═╡ 00ab445c-6c55-11eb-380f-2b64b223938e
md"# Isolating nodal line - proof

This notebook contains the code to prove that the nodal line for $u_2$ is isolated. It uses parameters for the curve $\Gamma$ computed in `4a-isolate-nodal-line-exploration.jl`.
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

# ╔═╡ 76bebd18-6ec1-11eb-2a8e-29e705c884ea
md"## Prove that $u_2$ is negative on $\Gamma$

Define a parameterization of the part of $\Gamma$ where we want to prove that $u_2$ is negative.
"

# ╔═╡ 57227a8e-6c56-11eb-3681-eb0a5b4e775d
p = let distance = domain.parent(load("../data/distance.jld")["distance"])
    y_max = distance * tanpi(domain.parent(1 // 6))
    p(t::Union{arb,arb_series}) = SVector{2,typeof(t)}(distance + 0t, t * y_max)
end

# ╔═╡ 51a735ee-6c57-11eb-1ff0-17b875452545
md"Compute an enclosure of the maximum value of $\tilde{u}_2(p(t))$ for $0 \leq t \leq 1$. Notice that $\tilde{u}_2$ is negative on $\Gamma$ so this will tell us how far away from zero it is."

# ╔═╡ a9ce8d12-6c57-11eb-0260-9d59be9e0d95
Γ_max = enclosemaximum(
	t -> u(p(t), λ),
	domain.parent(0),
	domain.parent(1),
	evaltype = :taylor,
	n = length(coefficients(u)) ÷ 8,
	rtol = 1e-2,
	lower_bound = u(p(domain.parent(0)), λ),
)

# ╔═╡ e09ab2fe-7601-11eb-3fc7-6b6afec9a531
md"The next step is to compute the $L^\infty$ bound for $\tilde{u}_2$. For this we need bounds for $\mu$, $\alpha$, $g(x)$ and $\max_{x \in \partial\Omega} |\tilde{u}_2(x)|$. For $g(x)$ we have the upper bound"

# ╔═╡ f0e495bc-6d46-11eb-34f1-278967a98791
g = sqrt(2MethodOfParticularSolutions.area(domain)) / 4domain.parent(π)

# ╔═╡ 14def084-6d47-11eb-2488-3d00b9ef8810
md"The lower bound for $\alpha$ we get by computing the distance between $\tilde{\lambda}_2$ and the interval $Λ'$ contain the third and fourth eigenvalue."

# ╔═╡ 4e61af5e-6d47-11eb-0a67-3d2b788d6aaf
α = let Λ´ = ArbTools.arb_load_dump(load("../data/cluster.jld")["Λ´_dump"], domain.parent)
    ArbTools.abs_lbound(Λ´ - λ)
end

# ╔═╡ bbfd810c-767d-11eb-1977-459b5bb31f9b
md"Both $\mu$ and the bound on the boundary have been computed in the notebook `3-isolating-second-jl` so we load those values"

# ╔═╡ d723e750-767d-11eb-0bc6-8793264e40d5
μ, m = let 
	μs_dump, max_boundary_dump = load("../data/enclosures.jld", "μs_dump", "max_boundary_dump")
	μ = ArbTools.arb_load_dump(μs_dump[2], domain.parent)
	m = ArbTools.arb_load_dump(max_boundary_dump[2], domain.parent)
	μ, m
end

# ╔═╡ bfa5c306-767e-11eb-1de3-1b955d68d8cd
md"This gives us the $L^\infty$ bound"

# ╔═╡ a073992a-767e-11eb-0faa-3188ad475f86
L_inf_bound = m * (1 + g * λ * (1 / (1 - μ) + 1 / α * (1 + μ^2 / α^2)))

# ╔═╡ c7c73892-767e-11eb-1f4f-8db3ab7b5022
md"If this bound is smaller than the absolute value of the maximum value on $\Gamma$, and this maximum value indeed is negative, then we can conclude that $u_2$ is negative on all of $\Gamma$."

# ╔═╡ dcd1e55c-767e-11eb-2164-4de6a1e7f8e0
u_is_negative_on_Γ = Γ_max < 0 && L_inf_bound < abs(Γ_max)

# ╔═╡ 8c1bf674-6ec3-11eb-2498-898f54c35a7b
md"## Prove that $u_2$ is positive on a points inside

For this we just have to evalute $\tilde{u}_2$ on a points inside and compare with the $L^\infty$ bound.
"

# ╔═╡ bb0b6fd2-6ec3-11eb-1b94-8f86f0b4710f
begin
    point_inside = [domain.parent(1 // 10), domain.parent(0)]
    u_on_point_inside = u(point_inside, λ)
    u_is_positive_inside = u_on_point_inside > L_inf_bound
end

# ╔═╡ 6e2ae6a6-6ec4-11eb-0bd8-9be009a1df4f
md"## Conclude

If $u_2$ is negative on $Γ$ and positive on a points inside $\Gamma$ then we can conclude that the nodal line is contained inside $\Gamma$."

# ╔═╡ 9afa04f0-6ec4-11eb-0bd7-f38df201dcaf
if u_is_negative_on_Γ && u_is_positive_inside
    md"Success! The nodal line is contained in $\Gamma$! 🎈🎈🎈"
else
    md"Could not conclude 😢"
end

# ╔═╡ 1ae695f6-6e9e-11eb-07dd-ff3d2f5aaecd
md"Finally we save the data that is used in the paper"

# ╔═╡ 2219744c-6e9e-11eb-0c30-99254c5ecf2a
let dir = "../data"
    save(
        joinpath(dir, "nodal-line.jld"),
        "Γ_max_dump",
        ArbTools.arb_dump(Γ_max),
        "α_dump",
        ArbTools.arb_dump(α),
        "L_infinity_bound_dump",
        ArbTools.arb_dump(L_inf_bound),
        "point_inside_dump",
        ArbTools.arb_dump.(point_inside),
        "u_on_point_inside_dump",
        ArbTools.arb_dump(u_on_point_inside),
    )
end

# ╔═╡ Cell order:
# ╟─00ab445c-6c55-11eb-380f-2b64b223938e
# ╠═5dbb9e30-6c55-11eb-36a3-210e1d52a1b5
# ╟─8039961a-6c55-11eb-2072-536d9ddf52db
# ╠═8784e15e-6c55-11eb-234c-eb1f0f423bae
# ╟─889a9cd0-6c55-11eb-2a2f-09e6225e11d4
# ╠═96c1da6e-6c55-11eb-3969-090f05e7647a
# ╟─76bebd18-6ec1-11eb-2a8e-29e705c884ea
# ╠═57227a8e-6c56-11eb-3681-eb0a5b4e775d
# ╟─51a735ee-6c57-11eb-1ff0-17b875452545
# ╠═a9ce8d12-6c57-11eb-0260-9d59be9e0d95
# ╟─e09ab2fe-7601-11eb-3fc7-6b6afec9a531
# ╠═f0e495bc-6d46-11eb-34f1-278967a98791
# ╟─14def084-6d47-11eb-2488-3d00b9ef8810
# ╠═4e61af5e-6d47-11eb-0a67-3d2b788d6aaf
# ╟─bbfd810c-767d-11eb-1977-459b5bb31f9b
# ╠═d723e750-767d-11eb-0bc6-8793264e40d5
# ╟─bfa5c306-767e-11eb-1de3-1b955d68d8cd
# ╠═a073992a-767e-11eb-0faa-3188ad475f86
# ╟─c7c73892-767e-11eb-1f4f-8db3ab7b5022
# ╠═dcd1e55c-767e-11eb-2164-4de6a1e7f8e0
# ╟─8c1bf674-6ec3-11eb-2498-898f54c35a7b
# ╠═bb0b6fd2-6ec3-11eb-1b94-8f86f0b4710f
# ╟─6e2ae6a6-6ec4-11eb-0bd8-9be009a1df4f
# ╠═9afa04f0-6ec4-11eb-0bd7-f38df201dcaf
# ╟─1ae695f6-6e9e-11eb-07dd-ff3d2f5aaecd
# ╠═2219744c-6e9e-11eb-0c30-99254c5ecf2a
