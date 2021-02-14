### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 8453cee0-6dd4-11eb-247a-a7ab8e8d7592
using ArbTools, MethodOfParticularSolutions, Nemo, PaynePolygon, Plots

# ╔═╡ 23d68114-6c76-11eb-373e-7b8358f16372
md"# Isolating the second eigenvalue

In this notebook we compute enclosures for the first to fourth eigenvalue, this then allows us to isolate the second eigenvalue.
"

# ╔═╡ b99db082-6dd4-11eb-2819-67e94af4dfe2
md"Load the precomputed approximations"

# ╔═╡ 4ef6f0ea-6dd4-11eb-094d-ebd29329adc8
domains, us, λs = begin
    domains_us_λs = [
        PaynePolygon.load_eigenfunction(
			"../data/approximate-eigenfunction-$i.jld", 
			T = arb,
		) 
		for i in 1:4
    ]

    tuple(zip(domains_us_λs...)...)
end

# ╔═╡ c4b6a6fe-6dd4-11eb-01ef-47b2c8121f1a
md"We need to compute enclosures of the eigenvalues for all four eigenfunctions. For this we need to bound them on the boundary and lower bound their norms. We first lower bound the norms"

# ╔═╡ 3bf6a08e-6dd5-11eb-056a-79ed2b8f6bc6
norms = [
	MethodOfParticularSolutions.norm(domains[i], us[i], λs[i]) for i in eachindex(domains)
]

# ╔═╡ 6147348e-6dd5-11eb-3cc3-6595d4e1fe57
md"For the boundary we don't need the tightest bound possible. Instead we compute an approximate bound by evaluating on a number of points and then we prove that they are bounded by a small multiple of this."

# ╔═╡ 8ca4b5c0-6dd5-11eb-105c-c3a165f91b93
approx_max = let
	compute_approx_max(domain, u, λ) = begin
		max_numpoints = length(coefficients(u))
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

# ╔═╡ 2be3d60c-6dd6-11eb-1510-bddeb7605885
Float64.(approx_max)

# ╔═╡ 5dffe32a-6e0c-11eb-0ac2-ff41ca816c6b
md"Not we prove that they are bounded by twice this approximate value on the boundary. For symmetry reasons we only have to bound them on a small subset of the boundary. **This is not quite what we do, we try to improve the bound until we get one which is smaller than twice the approximate value**."

# ╔═╡ 6882e56a-6e32-11eb-1117-efdda9116408
active_boundaries = [findall(!isempty, u.boundary_to_us) for u in us]

# ╔═╡ 46c2ca40-6e34-11eb-3560-5f084a97f76f
md"Below are the part of the boundary that we havea to check for the first/second and the third/fourth respectively."

# ╔═╡ 6768fefc-6e33-11eb-0e3e-6f7ead6a152c
let i = 1
	pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)
	for boundary in active_boundaries[i]
		p(t) = boundary_parameterization(t, domains[i], boundary)
		start = Float64.(p(0))
		stop = Float64.(p(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1)))
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
		stop = Float64.(p(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1)))
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

# ╔═╡ 64c58cfc-6e32-11eb-086c-2b90ab57e314


# ╔═╡ b315e674-6e0d-11eb-0c4d-49f4bac913d9
ok₁, res₁ = let i = 1
	ok = true
	res = domains[i].parent(0)
	for boundary in active_boundaries[i]
		stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1))
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

# ╔═╡ 8228d29a-6e35-11eb-3a2a-f3fa3c31e8a1
ok₂, res₂ = let i = 2
	ok = true
	res = domains[i].parent(0)
	for boundary in active_boundaries[i]
		stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1))
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

# ╔═╡ b48e4124-6e35-11eb-2321-bd2ab359d7c5
ok₃, res₃ = let i = 3
	ok = true
	res = domains[i].parent(0)
	for boundary in active_boundaries[i]
		stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1))
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
		stop = domains[i].parent(ifelse(boundary ∈ us[i].even_boundaries, 1//2, 1))
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
	md"Not all of the succeded!"
end

# ╔═╡ 0193502e-6e37-11eb-0959-0309ab0e8a81
max_values = [res₁, res₂, res₃, res₄]

# ╔═╡ 61ce8ef4-6e37-11eb-01fc-1f2c3988c964
md"Finally we can compute the enclosures!"

# ╔═╡ 662a78b4-6e37-11eb-376a-c1301458a078
enclosures = [
	begin
		ϵ = sqrt(MethodOfParticularSolutions.area(domains[i])) * max_values[i] / norms[i]	
		lower = λs[i] / (1 + getinterval(ϵ)[2])
    	upper = λs[i] / (1 - getinterval(ϵ)[2])
		setinterval(lower, upper)
	end
	for i in eachindex(domains)
]

# ╔═╡ 3f55f276-6e38-11eb-20a3-17477db84d0e
md"Check that the enclosures for the first and second eigenvalues are disjoint from the others."

# ╔═╡ 6a166b06-6e38-11eb-23ce-fb75d8ac8e71
all(!overlaps(enclosures[1], e) for e in enclosures[2:4])

# ╔═╡ b9e9a5d2-6e38-11eb-35c2-13de72c6d836
λ₂_isolated = all(!overlaps(enclosures[2], e) for e in [enclosures[1]; enclosures[3:4]])

# ╔═╡ ce629c12-6e38-11eb-2b9a-bb9abc6173d3
md"The third and fourth eigenvalue overlaps so instead we have to follow the procedure in the lemma in the paper."

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
md"For $g(x)$ we also have the bound"

# ╔═╡ d8672e3e-6e39-11eb-392f-8b6ec139f7da
g = sqrt(2MethodOfParticularSolutions.area(domains[1])) / 4domains[1].parent(π)

# ╔═╡ f1b47914-6e39-11eb-23ee-bd7dc6ef3387
md"Giving us"

# ╔═╡ f46d47ee-6e39-11eb-3369-610697b140a2
μ₃ = sqrt(MethodOfParticularSolutions.area(domains[3])) * max_values[3] / norms[3]

# ╔═╡ 19c3fe48-6e3a-11eb-2bb8-d56188827222
μ₄ = sqrt(MethodOfParticularSolutions.area(domains[4])) * max_values[4] / norms[4]

# ╔═╡ f95a5daa-6e39-11eb-18d7-b36d31b45ef4
md"and thus the $L^\infty$ bounds"

# ╔═╡ 2c9d1e64-6e3a-11eb-2bf6-71e8222ee966
bound₃ = max_values[3] * (1 + g * λs[3] * (1 / (1 - μ₃) + 1 / α * (1 + μ₃^2 / α^2)))

# ╔═╡ 48ee4a0c-6e3a-11eb-1d20-d9939ae5ab88
bound₄ = max_values[4] * (1 + g * λs[4] * (1 / (1 - μ₄) + 1 / α * (1 + μ₄^2 / α^2)))

# ╔═╡ 8b978df0-6e3a-11eb-075e-072925439fe5
md"We now evaluate both `u[3]` and `u[4]` at the points"

# ╔═╡ 951e7df4-6e3a-11eb-1c26-798f0b393333
points = [domains[3].parent.([0.5, 0.5]), domains[3].parent.([-0.5, 0.5])]

# ╔═╡ 46694af2-6e94-11eb-11b1-2bc77300a755
Float64[u(p, λ) for (u, λ) in ((us[3], λs[3]), (us[4], λs[4])), p in points]

# ╔═╡ 530d4704-6e3c-11eb-1103-c3cbe71e81d4
md"Check that we can definitely determine the signs at the points"

# ╔═╡ 668266e6-6e3c-11eb-2b65-699b9eb978f6
could_determine_sign = abs(us[3](points[1], λs[3])) - bound₃ > 0 &&
	abs(us[3](points[2], λs[3])) - bound₃ > 0 && 
	abs(us[4](points[1], λs[3])) - bound₄ > 0 &&
	abs(us[4](points[2], λs[3])) - bound₄ > 0

# ╔═╡ e514550a-6e3c-11eb-1724-439eb0841b50
md"Finally check that $u_3$ has different signs at the two points and that $u_4$ has the same sign"

# ╔═╡ f6c3fb8e-6e3c-11eb-1df7-95194b512b20
correct_signs = us[3](points[1], λs[3]) * us[3](points[2], λs[3]) < 0 &&
	us[4](points[1], λs[4]) * us[4](points[2], λs[4]) > 0

# ╔═╡ dd03a248-6e3d-11eb-20b3-71fdf50a05b0
md"To finish we check that everything actually succeded and so that the first and second eigenvalue do not overlap with $\Lambda'$"

# ╔═╡ cdb55222-6e93-11eb-1cbf-4b9f17d72c2b
begin
	λ₁_isolated = !overlaps(enclosures[1], enclosures[2]) && !overlaps(enclosures[1], Λ´)
	λ₁_isolated = !overlaps(enclosures[2], enclosures[1]) && !overlaps(enclosures[2], Λ´)
	if λ₁_isolated && λ₂_isolated && could_determine_sign && correct_signs
		md"Everything was successful! 🎈🎈🎈"
	else
		md"The proof did not succeed 😦"
	end
end

# ╔═╡ Cell order:
# ╟─23d68114-6c76-11eb-373e-7b8358f16372
# ╠═8453cee0-6dd4-11eb-247a-a7ab8e8d7592
# ╟─b99db082-6dd4-11eb-2819-67e94af4dfe2
# ╠═4ef6f0ea-6dd4-11eb-094d-ebd29329adc8
# ╟─c4b6a6fe-6dd4-11eb-01ef-47b2c8121f1a
# ╠═3bf6a08e-6dd5-11eb-056a-79ed2b8f6bc6
# ╟─6147348e-6dd5-11eb-3cc3-6595d4e1fe57
# ╠═8ca4b5c0-6dd5-11eb-105c-c3a165f91b93
# ╠═2be3d60c-6dd6-11eb-1510-bddeb7605885
# ╟─5dffe32a-6e0c-11eb-0ac2-ff41ca816c6b
# ╠═6882e56a-6e32-11eb-1117-efdda9116408
# ╟─46c2ca40-6e34-11eb-3560-5f084a97f76f
# ╠═6768fefc-6e33-11eb-0e3e-6f7ead6a152c
# ╟─26a68726-6e34-11eb-1d44-e5bf2700d36f
# ╟─64c58cfc-6e32-11eb-086c-2b90ab57e314
# ╠═b315e674-6e0d-11eb-0c4d-49f4bac913d9
# ╠═8228d29a-6e35-11eb-3a2a-f3fa3c31e8a1
# ╠═b48e4124-6e35-11eb-2321-bd2ab359d7c5
# ╠═c5d5a696-6e35-11eb-101c-4555a197b005
# ╟─15cb030c-6e37-11eb-2d26-c7676a668ed3
# ╠═0bd4cd1a-6e37-11eb-13a2-b11a437fb728
# ╠═0193502e-6e37-11eb-0959-0309ab0e8a81
# ╟─61ce8ef4-6e37-11eb-01fc-1f2c3988c964
# ╠═662a78b4-6e37-11eb-376a-c1301458a078
# ╟─3f55f276-6e38-11eb-20a3-17477db84d0e
# ╠═6a166b06-6e38-11eb-23ce-fb75d8ac8e71
# ╠═b9e9a5d2-6e38-11eb-35c2-13de72c6d836
# ╟─ce629c12-6e38-11eb-2b9a-bb9abc6173d3
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
# ╠═f46d47ee-6e39-11eb-3369-610697b140a2
# ╠═19c3fe48-6e3a-11eb-2bb8-d56188827222
# ╟─f95a5daa-6e39-11eb-18d7-b36d31b45ef4
# ╠═2c9d1e64-6e3a-11eb-2bf6-71e8222ee966
# ╠═48ee4a0c-6e3a-11eb-1d20-d9939ae5ab88
# ╟─8b978df0-6e3a-11eb-075e-072925439fe5
# ╠═951e7df4-6e3a-11eb-1c26-798f0b393333
# ╠═46694af2-6e94-11eb-11b1-2bc77300a755
# ╟─530d4704-6e3c-11eb-1103-c3cbe71e81d4
# ╠═668266e6-6e3c-11eb-2b65-699b9eb978f6
# ╟─e514550a-6e3c-11eb-1724-439eb0841b50
# ╠═f6c3fb8e-6e3c-11eb-1df7-95194b512b20
# ╟─dd03a248-6e3d-11eb-20b3-71fdf50a05b0
# ╠═cdb55222-6e93-11eb-1cbf-4b9f17d72c2b
