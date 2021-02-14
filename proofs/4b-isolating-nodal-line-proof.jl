### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 5dbb9e30-6c55-11eb-36a3-210e1d52a1b5
using ArbTools,
    JLD, MethodOfParticularSolutions, Nemo, NLsolve, PaynePolygon, Plots, StaticArrays

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

# ╔═╡ 76bebd18-6ec1-11eb-2a8e-29e705c884ea
md"## Prove that $u_2$ is negative on $\Gamma$

Define a parameterization of the part of $\Gamma$ where we want to prove that $u_2$ is negative.
"

# ╔═╡ 57227a8e-6c56-11eb-3681-eb0a5b4e775d
p = let distance = domain.parent(load("../data/distance.jld")["distance"])
    y_max = distance * tanpi(domain.parent(1 // 6))
    p(t::arb) = SVector{2,arb}(distance, t * y_max)
    # Hack to get the output of the correct type
    p(t::arb_series) = SVector{2,arb_series}(distance + 0t, t * y_max)
end

# ╔═╡ 51a735ee-6c57-11eb-1ff0-17b875452545
md"Compute an enclosure of the maximum value of `u(p(t), λ)` for `0 ≤ t ≤ 1`. Notice that `u` is negative on the line so this will tell us how far away from zero it is. We make use of the fact that the maximum in practice is attained at `t = 0` to speed up the computations."

# ╔═╡ a9ce8d12-6c57-11eb-0260-9d59be9e0d95
Γ_max = domain.parent(-1e-4) #enclosemaximum( # TODO
#t -> u(p(t), λ),
#domain.parent(0),
#domain.parent(1),
#evaltype = :taylor,
#n = length(coefficients(u)) ÷ 2,
#rtol = 1e-2,
#show_trace = true,
#lower_bound = u(p(domain.parent(0)), λ),
#)

# ╔═╡ 3427b152-6c58-11eb-073d-3b337dbf73a7
md"Once we have this value we can determine what $L^\infty$ bound we have to beat. Instead of computing the a as tight bound as possible we can then compute until we know the bound is good enough."

# ╔═╡ 7e4e1bd6-6c58-11eb-2ad1-958d24807648
md"We start by lower bounding $\alpha$, $g(x)$ and the norm of `u`."

# ╔═╡ f0e495bc-6d46-11eb-34f1-278967a98791
g = sqrt(2MethodOfParticularSolutions.area(domain)) / 4domain.parent(π)

# ╔═╡ 14def084-6d47-11eb-2488-3d00b9ef8810
md"The lower bound for $\alpha$ we get by computing the distance between `λ` and the interval `Λ´` contain the third and fourth eigenvalue."

# ╔═╡ 4e61af5e-6d47-11eb-0a67-3d2b788d6aaf
α = let Λ´ = ArbTools.arb_load_dump(load("../data/cluster.jld")["Λ´_dump"], domain.parent)
    ArbTools.abs_lbound(Λ´ - λ)
end

# ╔═╡ a94e6704-6d47-11eb-1b57-651063d065d2
n = MethodOfParticularSolutions.norm(domain, u, λ) # TODO: Fix this

# ╔═╡ 2c028666-6d49-11eb-0981-27a7f3e477c2
md"Now we can compute the $L^\infty$ bound given a bound `m` for $\max_{x \in \partial\Omega}|\tilde{u}_{2}(x)|$"

# ╔═╡ 588d6c28-6d49-11eb-21ed-bbdf23ddb0da
bound_from_max(m) = begin
    μ = sqrt(MethodOfParticularSolutions.area(domain)) * m / n
    return m * (1 + g * λ * (1 / (1 - μ) + 1 / α * (1 + μ^2 / α^2)))
end

# ╔═╡ 85944fac-6d49-11eb-3c04-4d4a49f20de7
md"We non-rigorously find `m` so that we get exactly `Γ_max`. We then take a value slightly lower by multiplying it with `0.99`, this is the bound that we'll prove that `u` satisfies."

# ╔═╡ 73868f1e-6d49-11eb-399c-3be79f711849
m_to_beat = let Γ_max = Float64(Γ_max)
    domain.parent(
        0.99nlsolve(m -> [Float64(bound_from_max(only(m))) - abs(Γ_max)], [Γ_max]).zero[1],
    )
end

# ╔═╡ 982fdaa4-6d4a-11eb-393b-514550552df2
md"Now rigorously check that we beat `Γ_max` if we can show that `u` is less than `m_to_beat`"

# ╔═╡ 461f3380-6d4b-11eb-3667-c5433ee84a8b
bound_from_max(m_to_beat) < Γ_max

# ╔═╡ 84a91dd0-6c58-11eb-224f-09806f358f5f
md"Now (assuming the above result is `true`) we are ready to prove that this bound indeed is satisfied on the boundary. For symmetry reasons there are only three boundaries for `domain` on which we have to bound `u`. We let `p₁`, `p₂` and `p₃` be parameterizations of them. Since `u` is symmetric on `p₁` and `p₃` we only have to bound the function on half of them."

# ╔═╡ ac61c2a8-6e94-11eb-0a2f-afd10d194562
begin
    active_boundaries = findall(!isempty, u.boundary_to_us)
    @assert length(active_boundaries) == 3
end

# ╔═╡ 28bf3b12-6d52-11eb-1e93-893daa638d64
let
    pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)
    for boundary in active_boundaries
        p(t) = boundary_parameterization(t, domain, boundary)
        start = Float64.(p(0))
        stop = Float64.(p(ifelse(boundary ∈ u.even_boundaries, 1 // 2, 1)))
        plot!(
            pl,
            [start[1], stop[1]],
            [start[2], stop[2]],
            color = :red,
            linewidth = 3,
            title = "Boundary to check",
        )
    end
    pl
end

# ╔═╡ d699905c-6d52-11eb-2739-d978535b5926
ok₁, res₁ = let boundary = active_boundaries[1]
    p(t) = boundary_parameterization(t, domain, boundary)
    stop = domain.parent(ifelse(boundary ∈ u.even_boundaries, 1 // 2, 1))
    PaynePolygon.bounded_by(
        t -> u(p(t), λ; boundary),
        domain.parent(stop), # TODO: Change this
        domain.parent(stop),
        m_to_beat,
        use_taylor = :true,
        n = length(coefficients(u)),
        start_intervals = 8,
        show_trace = true,
        show_evaluations = true,
        return_enclosure = true,
    )
end

# ╔═╡ 7506c67e-6d53-11eb-0533-832c62ac9929
ok₂, res₂ = let boundary = active_boundaries[2]
    p(t) = boundary_parameterization(t, domain, boundary)
    stop = domain.parent(ifelse(boundary ∈ u.even_boundaries, 1 // 2, 1))
    PaynePolygon.bounded_by(
        t -> u(p(t), λ; boundary),
        domain.parent(stop), # TODO: Change this 
        domain.parent(stop),
        m_to_beat,
        use_taylor = :true,
        n = length(coefficients(u)),
        start_intervals = 8,
        show_trace = true,
        show_evaluations = true,
        return_enclosure = true,
    )
end

# ╔═╡ 75b1bf02-6d53-11eb-360b-212da19f3faa
ok₃, res₃ = let boundary = active_boundaries[3]
    p(t) = boundary_parameterization(t, domain, boundary)
    stop = domain.parent(ifelse(boundary ∈ u.even_boundaries, 1 // 2, 1))
    PaynePolygon.bounded_by(
        t -> u(p(t), λ; boundary),
        domain.parent(stop), # TODO: Change this 
        domain.parent(stop),
        m_to_beat,
        use_taylor = :true,
        n = length(coefficients(u)),
        start_intervals = 8,
        show_trace = true,
        show_evaluations = true,
        return_enclosure = true,
    )
end

# ╔═╡ 9393479e-6dcf-11eb-3764-bd54c0d8393e
md"If all of these succeded then we have succesfully proved that $u_2$ is negative on all of $\Gamma$!"

# ╔═╡ 40027a24-6ec3-11eb-1df3-c5f1c5815dde
u_is_negative_on_Γ = ok₁ && ok₂ && ok₃

# ╔═╡ 2e6aedbc-6e98-11eb-041b-292ac90a3250
md"We also get that the bound on the boundary and the $L^\infty$ bound"

# ╔═╡ 38d73bfe-6e98-11eb-128e-bd9e3bcbb084
max_boundary = max(res₁, max(res₂, res₃))

# ╔═╡ ea55dbd2-6eba-11eb-1e9a-d79ea45ab0b5
L_infinity_bound = bound_from_max(max_boundary)

# ╔═╡ 8c1bf674-6ec3-11eb-2498-898f54c35a7b
md"## Prove that $u_2$ is positive on a points inside

For this we just have to evalute `u` on a points inside and compare with the $L^\infty$ bound.
"

# ╔═╡ bb0b6fd2-6ec3-11eb-1b94-8f86f0b4710f
begin
    point_inside = [domain.parent(1 // 10), domain.parent(0)]
    u_on_point_inside = u(point_inside, λ)
    u_is_positive_inside = u_on_point_inside > L_infinity_bound
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
        ArbTools.arb_dump(L_infinity_bound),
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
# ╟─3427b152-6c58-11eb-073d-3b337dbf73a7
# ╟─7e4e1bd6-6c58-11eb-2ad1-958d24807648
# ╠═f0e495bc-6d46-11eb-34f1-278967a98791
# ╟─14def084-6d47-11eb-2488-3d00b9ef8810
# ╠═4e61af5e-6d47-11eb-0a67-3d2b788d6aaf
# ╠═a94e6704-6d47-11eb-1b57-651063d065d2
# ╟─2c028666-6d49-11eb-0981-27a7f3e477c2
# ╠═588d6c28-6d49-11eb-21ed-bbdf23ddb0da
# ╟─85944fac-6d49-11eb-3c04-4d4a49f20de7
# ╠═73868f1e-6d49-11eb-399c-3be79f711849
# ╠═982fdaa4-6d4a-11eb-393b-514550552df2
# ╠═461f3380-6d4b-11eb-3667-c5433ee84a8b
# ╟─84a91dd0-6c58-11eb-224f-09806f358f5f
# ╠═ac61c2a8-6e94-11eb-0a2f-afd10d194562
# ╟─28bf3b12-6d52-11eb-1e93-893daa638d64
# ╠═d699905c-6d52-11eb-2739-d978535b5926
# ╠═7506c67e-6d53-11eb-0533-832c62ac9929
# ╠═75b1bf02-6d53-11eb-360b-212da19f3faa
# ╟─9393479e-6dcf-11eb-3764-bd54c0d8393e
# ╠═40027a24-6ec3-11eb-1df3-c5f1c5815dde
# ╟─2e6aedbc-6e98-11eb-041b-292ac90a3250
# ╠═38d73bfe-6e98-11eb-128e-bd9e3bcbb084
# ╠═ea55dbd2-6eba-11eb-1e9a-d79ea45ab0b5
# ╟─8c1bf674-6ec3-11eb-2498-898f54c35a7b
# ╠═bb0b6fd2-6ec3-11eb-1b94-8f86f0b4710f
# ╟─6e2ae6a6-6ec4-11eb-0bd8-9be009a1df4f
# ╠═9afa04f0-6ec4-11eb-0bd7-f38df201dcaf
# ╟─1ae695f6-6e9e-11eb-07dd-ff3d2f5aaecd
# ╠═2219744c-6e9e-11eb-0c30-99254c5ecf2a
