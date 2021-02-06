### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 55b3091a-67cf-11eb-0c5f-7d8286ce7fa7
using Arblib, LinearAlgebra, PaynePolygon

# ╔═╡ c44d75dc-67ce-11eb-21ae-ab7eebdcfc62
md"# Separating the first four eigenvalues

This contains the process for separating the first four eigenvalues, corresponding to section 3 in the paper.
"

# ╔═╡ 71b5010e-67cf-11eb-023e-f90a779af5ec
md"The domain is parametrised by three values `N`, `d` and `h`. The version in the article uses `N, d, h = 27, 11, 6` but to make this less computationally demanding we will use `9, 4, 3` as an example here."

# ╔═╡ a86a8ffe-67cf-11eb-0edb-e9ace2e00785
N, d, h = 9, 4, 3

# ╔═╡ cafaa066-67cf-11eb-23c4-136652f60d0e
md"As a first step we can plot the domain and the corresponding mesh."

# ╔═╡ dc02b16e-67cf-11eb-3546-73b6892586dd
PaynePolygon.plot_mesh(N, d, h)

# ╔═╡ 7007e2b2-67d0-11eb-14d7-955dc1932ff2
md"Around the holes the mesh looks slightly different."

# ╔═╡ 77b126a4-67d0-11eb-0340-47afc25c445f
PaynePolygon.plot_mesh(N, d, h, zoom_hole = true)

# ╔═╡ d17cfb9a-67d0-11eb-2325-35f7de366a98
md"The eigenvalue problem we want to solve is `Mx = λx` where `M = inv(B)A`. Here `A` is the stiffness matrix and `B` the mass matrix. The following code computes an enclosure of `M`."

# ╔═╡ 447e50c6-67d1-11eb-2e92-2168b2b7b36f
M = PaynePolygon.stiffness_matrix(Arb, N, d, h)

# ╔═╡ 8dde4d6e-67d1-11eb-237a-f5a30f822c55
md"We can compute approximations of the first five eigenvalues"

# ╔═╡ da5ee892-67dd-11eb-3e61-ed20b09df492
eigvals(convert(Hermitian{Float64,Matrix{Float64}}, M), 1:5)

# ╔═╡ 390b1e9c-67de-11eb-3fc5-a9733e66e0df
md"We now want to separate the first four eigenvalues from the rest using the method described in the paper. We first compute the `Q` matrix consisting of the approximate eigenvectors as columns."

# ╔═╡ c6817740-67ec-11eb-11a4-85eb3cf3ab41
Q = eigen(convert(Hermitian{Float64,Matrix{Float64}}, M)).vectors

# ╔═╡ 7bb6dfbe-6857-11eb-22d7-cd62481e93bc
md"The rest of the procedure is  implemented in `separate_eigenvalues` and when the `eltype` of `M` is `Arb` it returns rigorous results. We are guaranteed that `Λ` separates the first four eigenvalues of `M` from the rest."

# ╔═╡ 5c5eab8e-67de-11eb-1d7a-3b6e20d981d3
Λ = PaynePolygon.separate_eigenvalues(M, 4; Q)

# ╔═╡ Cell order:
# ╟─c44d75dc-67ce-11eb-21ae-ab7eebdcfc62
# ╠═55b3091a-67cf-11eb-0c5f-7d8286ce7fa7
# ╟─71b5010e-67cf-11eb-023e-f90a779af5ec
# ╠═a86a8ffe-67cf-11eb-0edb-e9ace2e00785
# ╟─cafaa066-67cf-11eb-23c4-136652f60d0e
# ╠═dc02b16e-67cf-11eb-3546-73b6892586dd
# ╟─7007e2b2-67d0-11eb-14d7-955dc1932ff2
# ╠═77b126a4-67d0-11eb-0340-47afc25c445f
# ╟─d17cfb9a-67d0-11eb-2325-35f7de366a98
# ╠═447e50c6-67d1-11eb-2e92-2168b2b7b36f
# ╟─8dde4d6e-67d1-11eb-237a-f5a30f822c55
# ╠═da5ee892-67dd-11eb-3e61-ed20b09df492
# ╟─390b1e9c-67de-11eb-3fc5-a9733e66e0df
# ╠═c6817740-67ec-11eb-11a4-85eb3cf3ab41
# ╟─7bb6dfbe-6857-11eb-22d7-cd62481e93bc
# ╠═5c5eab8e-67de-11eb-1d7a-3b6e20d981d3
