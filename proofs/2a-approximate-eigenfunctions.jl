### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 7ed715a0-6b7d-11eb-24c7-21588b2d267c
using Nemo, MethodOfParticularSolutions, PaynePolygon, Plots

# ╔═╡ 016dcf3e-685b-11eb-0b70-3f5dfe98193a
md"# Compute approximate eigenfunctions

This contains the process for computing the approximate eigenfunctions, corresponding to section 4 in the paper. The coefficients for the approximations are then stored and used in later computations.

We compute four eigenfunction $u_1$, $u_2$, $u_3$ and $u_4$, which require different levels of precision. All of them except $u_2$ can be computed to sufficient precision in a time that is reasonable for this type of notebook. For $u_2$ we compute it to lower precision than we in the end will need, the higher precision computations are done outside of this notebook.
"

# ╔═╡ a9c7d5e2-6b7d-11eb-1f70-996bf19341c9
import MethodOfParticularSolutions: example_domain_goal_v1

# ╔═╡ c2c6e4d2-6a26-11eb-2b0c-31de8f876374
md"We start by defining the four eigenfunctions we will use. Both $u_1$ and $u_2$ have a six-fold even symmetry, which we get by giving `symmetry_class = 1`. $u_3$ and $u_4$ have similar but slighly different symmetries, which we get with `symmetry_class = 2` and `symmetry_class = 3` respectively. The amount of precision required also differs between the eigenfunctions."

# ╔═╡ 345db9e6-6bac-11eb-1bad-7dd0c26b0449
begin
	symmetry_classes = [1, 1, 2, 3]
	precisions = [64, 384, 128, 128]
end

# ╔═╡ 659d4aa0-6b7d-11eb-12b5-f9adc7ba65cb
domains, us = let N = 27, d = 11, h = 6
    domains_us = [
        example_domain_goal_v1(
			N,
			d,
			h,
			ArbField(prec),
			T = Float64;
			symmetry_class,
			) for (prec, symmetry_class) in zip(precisions, symmetry_classes)
    ]

    tuple(zip(domains_us...)...)
end

# ╔═╡ 031b0290-6b7e-11eb-0924-8fae87c54ac2
md"We have the following approximations for the eigenvalues"

# ╔═╡ 0a9efd78-6b7e-11eb-1775-0db146884b8d
λs = [31.0432, 63.2083359862688427940, 63.7259, 63.7259]

# ╔═╡ 47263c60-6b7f-11eb-3d94-bb1a3a7f3690
md"The number of free coefficients we use is parametrised by `n`, for a given `n` the number of free coefficients used for the eigenfunction `u` will be `n*sum(u.orders)`. The values for `n` that we use are"

# ╔═╡ 5aefc90e-6b80-11eb-06c1-0ff8b3b4647d
ns = [1, 28, 6, 6]

# ╔═╡ ac3c2e10-6b80-11eb-2763-c3b7a7e3d582
md"And for references we have the following values for `sum(u.orders)`"

# ╔═╡ 9cdd5c42-6b81-11eb-2c36-6b248fc18a64
[sum(u.orders) for u in us]

# ╔═╡ fb3ee9ca-70f4-11eb-24ae-7fbacb3c583b
md"So the number of free coefficients for each eigenfunction is"

# ╔═╡ 0308bc62-70f5-11eb-21fe-7b14e9f01f7b
[n * sum(u.orders) for (n, u) in zip(ns, us)]

# ╔═╡ 1c997a18-7500-11eb-1b1f-555843dfa129
md"Another tuning parameter is the number of collocation points on the boundary (also the number of interior points used, but this rarely needs to be adjusted). The number of collocation points is given by a multiple of the number of free coefficients, this multiple is given below"

# ╔═╡ a65da562-7500-11eb-14cd-01a54e77c022
num_boundary_factor = [3, 16, 8, 8]

# ╔═╡ b964656e-6b81-11eb-05f2-79d153d02e39
md"We are now ready to compute the approximations. We compute the the approximation for each eigenfunction and then an approximate enclosure for the eigenvalue."

# ╔═╡ 07b535a0-6b80-11eb-37d2-5f39d69572db
let i = 1
	setprecision(BigFloat, precision(domains[i].parent)) do
    	mps!(
        	us[i],
        	domains[i],
        	BigFloat(λs[i]),
        	BigFloat(λs[i]),
        	ns[i] * sum(us[i].orders),
        	num_boundary = num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    	)
	end
end

# ╔═╡ f8cdf9c6-6b9d-11eb-3330-898302a734e5
let i = 1
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 4num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 954dabb8-6b99-11eb-25c9-550312317b24
let i = 2
	setprecision(BigFloat, precision(domains[i].parent)) do
    	mps!(
    	    us[i],
    	    domains[i],
    	    BigFloat(λs[i]),
    	    BigFloat(λs[i]),
    	    ns[i] * sum(us[i].orders),
    	    num_boundary = num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    	)
	end
end

# ╔═╡ 3917ba44-6b9e-11eb-33cb-55aba32d47c7
let i = 2
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 4num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 984bcbae-6b99-11eb-23fb-a33308ddfa83
let i = 3
	setprecision(BigFloat, precision(domains[i].parent)) do
    	mps!(
    	    us[i],
    	    domains[i],
    	    BigFloat(λs[i]),
    	    BigFloat(λs[i]),
    	    ns[i] * sum(us[i].orders),
    	    num_boundary = num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    	)
	end
end

# ╔═╡ 3be3d776-6b9e-11eb-3ba1-79a7b11659f0
let i = 3
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 4num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 9b0f97a0-6b99-11eb-1fdc-31553d6a2bbc
let i = 4
	setprecision(BigFloat, precision(domains[i].parent)) do
    	mps!(
       		us[i],
        	domains[i],
        	BigFloat(λs[i]),
        	BigFloat(λs[i]),
        	ns[i] * sum(us[i].orders),
        	num_boundary = num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    	)
	end
end

# ╔═╡ 3fc77686-6b9e-11eb-0ad9-0d7f0000e326
let i = 4
    MethodOfParticularSolutions.enclose_eigenvalue_approx(
        domains[i],
        us[i],
        domains[i].parent(λs[i]),
        max_numpoints = 4num_boundary_factor[i]*ns[i] * sum(us[i].orders),
    )
end

# ╔═╡ 17a0a64e-6b87-11eb-38e6-f755f2657de8
md"To be able to reuse the results later we want to store the computed coefficients. We can use `PaynePolygon.save_eigenfunction` for saving the required data to a JLD-file."

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
```julia
domain, u, λ = PaynePolygon.load_eigenfunction(\"../data/approximate-eigenfunction-$i.jld\")
```
"

# ╔═╡ ee251ab0-6b98-11eb-014c-b5eb069d786e
md"See `proofs/2b-approximate-eigenfunctions-plots.jl` for plots of the eigenfunctions."

# ╔═╡ Cell order:
# ╟─016dcf3e-685b-11eb-0b70-3f5dfe98193a
# ╠═7ed715a0-6b7d-11eb-24c7-21588b2d267c
# ╠═a9c7d5e2-6b7d-11eb-1f70-996bf19341c9
# ╟─c2c6e4d2-6a26-11eb-2b0c-31de8f876374
# ╠═345db9e6-6bac-11eb-1bad-7dd0c26b0449
# ╠═659d4aa0-6b7d-11eb-12b5-f9adc7ba65cb
# ╟─031b0290-6b7e-11eb-0924-8fae87c54ac2
# ╠═0a9efd78-6b7e-11eb-1775-0db146884b8d
# ╟─47263c60-6b7f-11eb-3d94-bb1a3a7f3690
# ╠═5aefc90e-6b80-11eb-06c1-0ff8b3b4647d
# ╟─ac3c2e10-6b80-11eb-2763-c3b7a7e3d582
# ╠═9cdd5c42-6b81-11eb-2c36-6b248fc18a64
# ╟─fb3ee9ca-70f4-11eb-24ae-7fbacb3c583b
# ╠═0308bc62-70f5-11eb-21fe-7b14e9f01f7b
# ╟─1c997a18-7500-11eb-1b1f-555843dfa129
# ╠═a65da562-7500-11eb-14cd-01a54e77c022
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
