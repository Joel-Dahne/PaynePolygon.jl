### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 7ec30e64-7061-11eb-30c0-11588d34f045
using Revise

# ╔═╡ e800aea2-705e-11eb-37d1-c51771f6f70b
using MethodOfParticularSolutions, PaynePolygon, Plots

# ╔═╡ 942342c0-705e-11eb-28c6-b17eb6e86323
md"# Plots of approximate eigenfunctions

The approximate eigenfunctions were computed in the previous notebook and we here give some plots of the approximations together with some other related plots.
"

# ╔═╡ 11869a04-7062-11eb-1b79-dd45bea5b782
import MethodOfParticularSolutions: vertices, center

# ╔═╡ f2626688-705e-11eb-02ab-0520ae1cbc31
domains, us, λs = let
    domains_us_λs = [
        PaynePolygon.load_eigenfunction(
            "../data/approximate-eigenfunction-$i.jld",
            T = Float64,
        ) for i = 1:4
    ]

    tuple(zip(domains_us_λs...)...)
end

# ╔═╡ 5a3367a6-7061-11eb-3c48-b937c8f9b006
let domain = domains[1]
	pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)
	
	scatter!(
		pl,
		Float64.(getindex.(vertices(domain.exterior), 1)),
		Float64.(getindex.(vertices(domain.exterior), 2)),
		color = :red,
		markersize = 7,
	)
	
	vs = [vertex(d, 2) for d in domain.interiors]
	scatter!(
		pl,
		Float64.(getindex.(vs, 1)),
		Float64.(getindex.(vs, 2)),
		color = :blue,
		markersize = 7,
	)
	
	vs = [vertex(d, i) for d in domain.interiors, i in [1, 3]][:]
	scatter!(
		pl,
		Float64.(getindex.(vs, 1)),
		Float64.(getindex.(vs, 2)),
		color = :green,
		markersize = 7,
	)
	
	scatter!(
		pl,
		Float64[MethodOfParticularSolutions.center(domain)[1]],
		Float64[MethodOfParticularSolutions.center(domain)[2]],
		color = :yellow,
		markersize = 7,
	)
end

# ╔═╡ 767111ce-706a-11eb-3245-49883691e3d8
N = 200

# ╔═╡ 1ace661c-7064-11eb-2b08-c92eafc00779
pl1 = let i = 1
	domain, u, λ = domains[i], us[i], λs[i]
	
	PaynePolygon.plot_eigenfunction(domain, u, λ, N, N)
end

# ╔═╡ 764ff670-7066-11eb-259f-8b6a17f2f122
pl2 = let i = 2
	domain, u, λ = domains[i], us[i], λs[i]
	
	PaynePolygon.plot_eigenfunction(domain, u, λ, N, N)
end

# ╔═╡ 0f902fe8-7068-11eb-3a97-914409e9fdf3
pl3 = let i = 3
	domain, u, λ = domains[i], us[i], λs[i]
	
	PaynePolygon.plot_eigenfunction(domain, u, λ, N, N)
end

# ╔═╡ 8147efce-706b-11eb-3ca2-6b4edc35b3a4
pl4 = let i = 4
	domain, u, λ = domains[i], us[i], λs[i]
	
	PaynePolygon.plot_eigenfunction(domain, u, λ, N, N)
end

# ╔═╡ d354203a-706b-11eb-3378-5f6fe7f05596
pl_all = plot(pl1, pl2, pl3, pl4)

# ╔═╡ b01014d6-706a-11eb-1bd4-95c246347d8a
let dir = "../figures"
	for (i, pl) in enumerate([pl1, pl2, pl3, pl4])
		savefig(pl, joinpath(dir, "approximate-eigenfunction-$i.pdf"))
	end
	savefig(pl_all, joinpath(dir, "approximate-eigenfunctions.pdf"))
end

# ╔═╡ Cell order:
# ╟─942342c0-705e-11eb-28c6-b17eb6e86323
# ╠═7ec30e64-7061-11eb-30c0-11588d34f045
# ╠═e800aea2-705e-11eb-37d1-c51771f6f70b
# ╠═11869a04-7062-11eb-1b79-dd45bea5b782
# ╠═f2626688-705e-11eb-02ab-0520ae1cbc31
# ╠═5a3367a6-7061-11eb-3c48-b937c8f9b006
# ╟─767111ce-706a-11eb-3245-49883691e3d8
# ╟─1ace661c-7064-11eb-2b08-c92eafc00779
# ╟─764ff670-7066-11eb-259f-8b6a17f2f122
# ╟─0f902fe8-7068-11eb-3a97-914409e9fdf3
# ╟─8147efce-706b-11eb-3ca2-6b4edc35b3a4
# ╟─d354203a-706b-11eb-3378-5f6fe7f05596
# ╠═b01014d6-706a-11eb-1bd4-95c246347d8a
