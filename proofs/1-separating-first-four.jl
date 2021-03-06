### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 55b3091a-67cf-11eb-0c5f-7d8286ce7fa7
using Arblib, GenericLinearAlgebra, JLD, LinearAlgebra, MultiFloats, PaynePolygon, Plots

# ╔═╡ c44d75dc-67ce-11eb-21ae-ab7eebdcfc62
md"# Separating the first four eigenvalues

This contains the process for separating the first four eigenvalues, corresponding to section 3 in the paper.
"

# ╔═╡ 8ce329d4-b48b-11eb-1c7d-1f6c3e14e91a
md"Set the precision used in the computations."

# ╔═╡ 42bb750a-b48b-11eb-270a-8937d22533cf
begin
    setprecision(BigFloat, 128)
    setprecision(Arb, 128)
end

# ╔═╡ 71b5010e-67cf-11eb-023e-f90a779af5ec
md"The domain is parametrised by three values `N`, `d` and `h`. The version in the article uses `N, d, h = 27, 11, 6`. However these computation require up to 32GB of memory so to make it less memory heavy you can use `9, 4, 3` for a demo."

# ╔═╡ 0fab0e74-74fc-11eb-3257-d75a61de832e
demo = true

# ╔═╡ a86a8ffe-67cf-11eb-0edb-e9ace2e00785
N, d, h = ifelse(demo, (9, 4, 3), (27, 11, 6))

# ╔═╡ cafaa066-67cf-11eb-23c4-136652f60d0e
md"As a first step we can plot the domain both with and without its mesh"

# ╔═╡ dc02b16e-67cf-11eb-3546-73b6892586dd
pl_domain = PaynePolygon.plot_mesh(N, d, h, plot_mesh = false)

# ╔═╡ d636096e-6c81-11eb-118f-797bdad6d403
pl_mesh = PaynePolygon.plot_mesh(N, d, h)

# ╔═╡ 7007e2b2-67d0-11eb-14d7-955dc1932ff2
md"Around the holes the mesh looks slightly different."

# ╔═╡ 77b126a4-67d0-11eb-0340-47afc25c445f
pl_mesh_zoomed = PaynePolygon.plot_mesh(N, d, h, zoom_hole = true)

# ╔═╡ f4785650-6c81-11eb-3845-d7d63013ca2b
md"We save these three figures for inclusion in the paper"

# ╔═╡ 019eb77c-6c82-11eb-1f06-cf9705b47624
let dir = "../figures"
    savefig(pl_domain, joinpath(dir, "domain.pdf"))
    savefig(pl_mesh, joinpath(dir, "mesh.pdf"))
    savefig(pl_mesh_zoomed, joinpath(dir, "mesh-zoomed.pdf"))
end

# ╔═╡ d17cfb9a-67d0-11eb-2325-35f7de366a98
md"The eigenvalue problem we want to solve is $Mx = λx$ where $M = B^{-1/2}AB^{-1/2}$. Here $A$ is the stiffness matrix and $B$ the mass matrix. Since $A$ is symmetric and $B$ is diagonal we have that $M$ is symmetric. The following code computes an enclosure of $M$."

# ╔═╡ 447e50c6-67d1-11eb-2e92-2168b2b7b36f
M = PaynePolygon.stiffness_matrix(Arb, N, d, h)

# ╔═╡ 8dde4d6e-67d1-11eb-237a-f5a30f822c55
md"We can compute approximations of the first five eigenvalues"

# ╔═╡ da5ee892-67dd-11eb-3e61-ed20b09df492
eigvals(LinearAlgebra.copy_oftype(M, Float64), 1:5)

# ╔═╡ 390b1e9c-67de-11eb-3fc5-a9733e66e0df
md"We now want to separate the first four eigenvalues from the rest using the method described in the paper. We first compute the $Q$ matrix consisting of the approximate eigenvectors as columns."

# ╔═╡ dd74cc60-6c84-11eb-23d9-cb99aecdff47
md"Unfortunately we need the eigenvectors to slightly higher that `Float64` precision for the separation to succed in the case when `N, d, h = 27, 11, 6`. To accomplish that we use `Float64x2` from [MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl) which gives us twice the precision. We then also have to use a generic version of `eigen` because the one implemented in `LinearAlgebra` only supports machine floats, this is implemented in [GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl)."

# ╔═╡ bf2a356e-6c8a-11eb-25ee-97dccf30de05
md"We need two helper methods for converting between `Arb` and `Float64x2`"

# ╔═╡ a35df8d4-6c87-11eb-0fb6-715d86edc1a7
begin
    MultiFloats.MultiFloat{Float64,N}(x::Arb) where {N} =
        MultiFloat{Float64,N}(BigFloat(Arblib.midref(x)))
    Arblib.set!(res::Arblib.ArbLike, x::MultiFloat{Float64,N}) where {N} =
        Arblib.set!(res, BigFloat(x))
end

# ╔═╡ e5805844-6e06-11eb-1739-ff95933be386
md"Since this is one of the more costly parts of the computations we have worked to speed up the `eigen` version from GenericLinearAlgebra.jl. The most costly part of the computation is conversion to a symmetric tridiagonal matrix, this part we have parallelized and improved slightly in `PaynePolygon.symtri!`, we therefore manually perform this operation first."

# ╔═╡ 8737620e-6e07-11eb-3888-7de505b82886
Q = let
    A = GenericLinearAlgebra.SymmetricTridiagonalFactorization(
        PaynePolygon.symtri!(LinearAlgebra.copy_oftype(M, Float64x2))...,
    )
    B = PaynePolygon._Array(A.Q)

    ArbMatrix(
        LinearAlgebra.Eigen(
            GenericLinearAlgebra.eigQL!(
                A.diagonals,
                vectors = B,
                tol = convert(eltype(B), 1e-20),
            )...,
        ).vectors,
    )
end;

# ╔═╡ 7bb6dfbe-6857-11eb-22d7-cd62481e93bc
md"The rest of the procedure is  implemented in `separate_eigenvalues` and when the `eltype` of `M` is `Arb` it returns rigorous results. It returns `s, λ_h_4_upper, λ_h_5_lower` where `s` is the $s$ in Lemma 3.2 in the paper, `λ_h_4_upper` is an upper bound of the fourth eigenvalue and `λ_h_5_lower` is a lower bound of the fifth eigenvalue. We are guaranteed that `λ_h_5_lower` separates the first four eigenvalues of `M` from the rest."

# ╔═╡ 5c5eab8e-67de-11eb-1d7a-3b6e20d981d3
s, λ_h_4_upper, λ_h_5_lower = PaynePolygon.separate_eigenvalues(M, 4; Q)

# ╔═╡ 9a794330-712e-11eb-34f9-5fd98af4fcbe
md"We can now use Theorem 3.1 from the paper to get a lower bound of the fifth eigenvalue for the eigenvalue problem"

# ╔═╡ b80114de-712e-11eb-33c9-4bb0e7f354e8
λ_5_lower = let C = Arb("0.1893") / N
    λ_h_5_lower / (1 + C^2 * λ_h_5_lower)
end

# ╔═╡ 10f7f10a-6a09-11eb-16e7-d5379bcfd230
md"Finally we save the results"

# ╔═╡ ab9ac1f6-6c82-11eb-3ea0-7574d2ada4d2
save(
    "../data/separation-bound.jld",
    "λ_h_5_lower_dump",
    Arblib.dump_string(λ_h_5_lower),
    "λ_5_lower_dump",
    Arblib.dump_string(λ_5_lower),
)

# ╔═╡ Cell order:
# ╟─c44d75dc-67ce-11eb-21ae-ab7eebdcfc62
# ╠═55b3091a-67cf-11eb-0c5f-7d8286ce7fa7
# ╟─8ce329d4-b48b-11eb-1c7d-1f6c3e14e91a
# ╠═42bb750a-b48b-11eb-270a-8937d22533cf
# ╟─71b5010e-67cf-11eb-023e-f90a779af5ec
# ╠═0fab0e74-74fc-11eb-3257-d75a61de832e
# ╠═a86a8ffe-67cf-11eb-0edb-e9ace2e00785
# ╟─cafaa066-67cf-11eb-23c4-136652f60d0e
# ╟─dc02b16e-67cf-11eb-3546-73b6892586dd
# ╟─d636096e-6c81-11eb-118f-797bdad6d403
# ╟─7007e2b2-67d0-11eb-14d7-955dc1932ff2
# ╟─77b126a4-67d0-11eb-0340-47afc25c445f
# ╟─f4785650-6c81-11eb-3845-d7d63013ca2b
# ╠═019eb77c-6c82-11eb-1f06-cf9705b47624
# ╟─d17cfb9a-67d0-11eb-2325-35f7de366a98
# ╠═447e50c6-67d1-11eb-2e92-2168b2b7b36f
# ╟─8dde4d6e-67d1-11eb-237a-f5a30f822c55
# ╠═da5ee892-67dd-11eb-3e61-ed20b09df492
# ╟─390b1e9c-67de-11eb-3fc5-a9733e66e0df
# ╟─dd74cc60-6c84-11eb-23d9-cb99aecdff47
# ╟─bf2a356e-6c8a-11eb-25ee-97dccf30de05
# ╠═a35df8d4-6c87-11eb-0fb6-715d86edc1a7
# ╟─e5805844-6e06-11eb-1739-ff95933be386
# ╠═8737620e-6e07-11eb-3888-7de505b82886
# ╟─7bb6dfbe-6857-11eb-22d7-cd62481e93bc
# ╠═5c5eab8e-67de-11eb-1d7a-3b6e20d981d3
# ╟─9a794330-712e-11eb-34f9-5fd98af4fcbe
# ╠═b80114de-712e-11eb-33c9-4bb0e7f354e8
# ╟─10f7f10a-6a09-11eb-16e7-d5379bcfd230
# ╠═ab9ac1f6-6c82-11eb-3ea0-7574d2ada4d2
