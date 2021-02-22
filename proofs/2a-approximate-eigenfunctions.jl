### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 7ed715a0-6b7d-11eb-24c7-21588b2d267c
using MethodOfParticularSolutions, PaynePolygon, Plots

# ╔═╡ 016dcf3e-685b-11eb-0b70-3f5dfe98193a
md"# Computing approximate eigenfunctions

This contains the process for computing the approximate eigenfunctions, corresponding to section 4 in the paper. The coefficients for the approximations are then stored and used in later computations.

We compute four eigenfunction `u₁`, `u₂`, `u₃` and `u₄`, which require different levels of precision. All of them except `u₂` can be computed to sufficient precision in a time that is reasonable for this type of notebook. For `u₂` we compute it to lower precision than we in the end will need, the higher precision computations are done outside of this notebook.
"

# ╔═╡ a9c7d5e2-6b7d-11eb-1f70-996bf19341c9
import MethodOfParticularSolutions: example_domain_goal_v1

# ╔═╡ c2c6e4d2-6a26-11eb-2b0c-31de8f876374
md"We start by defining the four eigenfunctions we will use. Both $u_1$ and $u_2$ have a six-fold even symmetry, which we get by giving `symmetry_class = 1`. $u_3$ and $u_4$ have similar but slighly different symmetries, which we get with `symmetry_class = 2` and `symmetry_class = 3` respectively."

# ╔═╡ 345db9e6-6bac-11eb-1bad-7dd0c26b0449
symmetry_classes = [1, 1, 2, 3]

# ╔═╡ 659d4aa0-6b7d-11eb-12b5-f9adc7ba65cb
domains, us = let N = 27, d = 11, h = 6
    domains_us = [
        example_domain_goal_v1(N, d, h, T = Float64; symmetry_class) for
        symmetry_class in symmetry_classes
    ]

    tuple(zip(domains_us...)...)
end

# ╔═╡ 031b0290-6b7e-11eb-0924-8fae87c54ac2
md"We have the following approximations for the eigenvalues"

# ╔═╡ 0a9efd78-6b7e-11eb-1775-0db146884b8d
λs = [31.0432, 63.20833838, 63.7259, 63.7259]

# ╔═╡ 1701c14e-6b7f-11eb-27ef-dd57dff8f572
md"We can plot the domain that is being used"

# ╔═╡ 1de365da-6b7f-11eb-191b-ed48b4638e1a
plot(domains[1], 1000, 0, legend = :none)

# ╔═╡ 47263c60-6b7f-11eb-3d94-bb1a3a7f3690
md"Finally we want to compute the approximations. The number of free coefficients we use is parametrised by `n`, for a given `n` the number of free coefficients used for the eigenfunction `u` will be `n*sum(u.orders)`. The values for `n` that we use are"

# ╔═╡ 5aefc90e-6b80-11eb-06c1-0ff8b3b4647d
ns = [1, 10, 6, 6]

# ╔═╡ ac3c2e10-6b80-11eb-2763-c3b7a7e3d582
md"And for references we have the following values for `sum(u.orders)`"

# ╔═╡ 9cdd5c42-6b81-11eb-2c36-6b248fc18a64
[sum(u.orders) for u in us]

# ╔═╡ fb3ee9ca-70f4-11eb-24ae-7fbacb3c583b
md"So the number of free coefficients for each eigenfunction is"

# ╔═╡ 0308bc62-70f5-11eb-21fe-7b14e9f01f7b
[n * sum(u.orders) for (n, u) in zip(ns, us)]

# ╔═╡ b964656e-6b81-11eb-05f2-79d153d02e39
md"We are now ready to compute the approximations. We compute the the approximation for each eigenfunction and then an approximate enclosure for the eigenvaleu."

# ╔═╡ 07b535a0-6b80-11eb-37d2-5f39d69572db
let i = 1
    mps!(
        us[i],
        domains[i],
        BigFloat(λs[i]),
        BigFloat(λs[i]),
        ns[i] * sum(us[i].orders),
        num_boundary = 3ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ f8cdf9c6-6b9d-11eb-3330-898302a734e5
let i = 1
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 3ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 954dabb8-6b99-11eb-25c9-550312317b24
let i = 2
    mps!(
        us[i],
        domains[i],
        BigFloat(λs[i]),
        BigFloat(λs[i]),
        ns[i] * sum(us[i].orders),
        num_boundary = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 3917ba44-6b9e-11eb-33cb-55aba32d47c7
let i = 2
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 984bcbae-6b99-11eb-23fb-a33308ddfa83
let i = 3
    mps!(
        us[i],
        domains[i],
        BigFloat(λs[i]),
        BigFloat(λs[i]),
        ns[i] * sum(us[i].orders),
        num_boundary = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 3be3d776-6b9e-11eb-3ba1-79a7b11659f0
let i = 3
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 9b0f97a0-6b99-11eb-1fdc-31553d6a2bbc
let i = 4
    mps!(
        us[i],
        domains[i],
        BigFloat(λs[i]),
        BigFloat(λs[i]),
        ns[i] * sum(us[i].orders),
        num_boundary = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 3fc77686-6b9e-11eb-0ad9-0d7f0000e326
let i = 4
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 8ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 17a0a64e-6b87-11eb-38e6-f755f2657de8
md"To be able to reuse the results later we want to store the computed coefficients. We can use `PaynePolygon.save_eigenfunctions` for saving the required data to a JLD-file."

# ╔═╡ 2d2342d6-6b98-11eb-2f60-a9b357272958
for i in eachindex(us)
    PaynePolygon.save_eigenfunction(
        "../data/approximate-eigenfunction-$i.jld",
        us[i],
        λs[i],
        symmetry_class = symmetry_classes[i],
    )
end

# ╔═╡ e5e2103e-6b88-11eb-3aff-61834843c2ad
md"The eigenfunctions can then be reconstructed at a later time with
```
domain, u, λ = PaynePolygon.load_eigenfunctions(\"../data/approximate-eigenfunctions-$i.jld\")
```
"

# ╔═╡ ee251ab0-6b98-11eb-014c-b5eb069d786e
md"We end by plotting the eigenfunctions"

# ╔═╡ fce2f214-6b98-11eb-075d-976608a5f0fb
eigenfunctionheatmap(domains[1], us[1], λs[1], 50, 50)

# ╔═╡ 30b275f6-6b99-11eb-30b2-4f440492de8c
eigenfunctionheatmap(domains[2], us[2], λs[2], 50, 50)

# ╔═╡ 3d713886-6b99-11eb-0eb6-05c498c8ea2b
eigenfunctionheatmap(domains[3], us[3], λs[3], 50, 50)

# ╔═╡ 4146eaa0-6b99-11eb-2411-67e6d5cf30fe
eigenfunctionheatmap(domains[4], us[4], λs[4], 50, 50)

# ╔═╡ 30dc3c2a-6b9d-11eb-1691-2f73520e669b
md"We can highlight the nodal line of the second eigenfunction. **TODO**: Make a better plot of this."

# ╔═╡ 8159d7ca-6b9d-11eb-0be5-b321ee4c15f1
eigenfunctionheatmap(
    domains[2],
    us[2],
    λs[2],
    range(-0.5, 0.5, length = 100),
    range(-0.5, 0.5, length = 100),
    absolute_value = true,
)

# ╔═╡ Cell order:
# ╟─016dcf3e-685b-11eb-0b70-3f5dfe98193a
# ╠═7ed715a0-6b7d-11eb-24c7-21588b2d267c
# ╠═a9c7d5e2-6b7d-11eb-1f70-996bf19341c9
# ╟─c2c6e4d2-6a26-11eb-2b0c-31de8f876374
# ╠═345db9e6-6bac-11eb-1bad-7dd0c26b0449
# ╠═659d4aa0-6b7d-11eb-12b5-f9adc7ba65cb
# ╟─031b0290-6b7e-11eb-0924-8fae87c54ac2
# ╠═0a9efd78-6b7e-11eb-1775-0db146884b8d
# ╟─1701c14e-6b7f-11eb-27ef-dd57dff8f572
# ╠═1de365da-6b7f-11eb-191b-ed48b4638e1a
# ╟─47263c60-6b7f-11eb-3d94-bb1a3a7f3690
# ╠═5aefc90e-6b80-11eb-06c1-0ff8b3b4647d
# ╟─ac3c2e10-6b80-11eb-2763-c3b7a7e3d582
# ╠═9cdd5c42-6b81-11eb-2c36-6b248fc18a64
# ╠═fb3ee9ca-70f4-11eb-24ae-7fbacb3c583b
# ╠═0308bc62-70f5-11eb-21fe-7b14e9f01f7b
# ╟─b964656e-6b81-11eb-05f2-79d153d02e39
# ╠═07b535a0-6b80-11eb-37d2-5f39d69572db
# ╠═f8cdf9c6-6b9d-11eb-3330-898302a734e5
# ╠═954dabb8-6b99-11eb-25c9-550312317b24
# ╠═3917ba44-6b9e-11eb-33cb-55aba32d47c7
# ╠═984bcbae-6b99-11eb-23fb-a33308ddfa83
# ╠═3be3d776-6b9e-11eb-3ba1-79a7b11659f0
# ╠═9b0f97a0-6b99-11eb-1fdc-31553d6a2bbc
# ╠═3fc77686-6b9e-11eb-0ad9-0d7f0000e326
# ╟─17a0a64e-6b87-11eb-38e6-f755f2657de8
# ╠═2d2342d6-6b98-11eb-2f60-a9b357272958
# ╟─e5e2103e-6b88-11eb-3aff-61834843c2ad
# ╟─ee251ab0-6b98-11eb-014c-b5eb069d786e
# ╠═fce2f214-6b98-11eb-075d-976608a5f0fb
# ╠═30b275f6-6b99-11eb-30b2-4f440492de8c
# ╠═3d713886-6b99-11eb-0eb6-05c498c8ea2b
# ╠═4146eaa0-6b99-11eb-2411-67e6d5cf30fe
# ╟─30dc3c2a-6b9d-11eb-1691-2f73520e669b
# ╠═8159d7ca-6b9d-11eb-0be5-b321ee4c15f1
