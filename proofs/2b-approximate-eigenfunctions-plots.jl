### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

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

# ╔═╡ bbb33f6a-7769-11eb-1275-b74033697b16
md"First we give a plot showing the placement of the basis functions for the computed approximations. For the first and second eigenfunction the points with the same colors correspond to basis functions with the same coefficients. For the third and fourth eigenfunction the symmetries are however different. We save this to `../figures/placement.pdf`."

# ╔═╡ 5a3367a6-7061-11eb-3c48-b937c8f9b006
let domain = domains[1]
    pl = PaynePolygon.plot_mesh(27, 11, 6, plot_mesh = false)

    scatter!(
        pl,
        Float64.(getindex.(vertices(domain.exterior), 1)),
        Float64.(getindex.(vertices(domain.exterior), 2)),
		label = "First type",
		legendfontsize = 11,
        color = :red,
        markersize = 7,
    )

    vs = [vertex(d, 2) for d in domain.interiors]
    scatter!(
        pl,
        Float64.(getindex.(vs, 1)),
        Float64.(getindex.(vs, 2)),
		label = "Second type",
        color = :blue,
        markersize = 7,
    )

    vs = [vertex(d, i) for d in domain.interiors, i in [1, 3]][:]
    scatter!(
        pl,
        Float64.(getindex.(vs, 1)),
        Float64.(getindex.(vs, 2)),
		label = "Second type",
        color = :green,
        markersize = 7,
    )

    scatter!(
        pl,
        Float64[MethodOfParticularSolutions.center(domain)[1]],
        Float64[MethodOfParticularSolutions.center(domain)[2]],
		label = "Third type",
        color = :yellow,
        markersize = 7,
    )
	
	# Charges
	for u in [us[1].us[2].us; us[1].us[3].us]
		charges = let n = 5
        	cs = [MethodOfParticularSolutions.charge(u, i, n, true) for i = 1:n]
        	(Float64.(getindex.(cs, 1)), Float64.(getindex.(cs, 2)))
    	end
		scatter!(charges[1], charges[2], label = "", color = :red, markersize = 2)
	end
	
	savefig(pl, "../figures/placement.pdf")
	
	pl
end

# ╔═╡ e1c8c5a2-776a-11eb-3a87-97f694a5429a
md"Next we simply plot the four approximate eigenfunctions."

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

# ╔═╡ f51966c8-776a-11eb-2667-f5671d29a85f
md"Finally we make a plot higlighting the nodal line of the second eigenfunction. This is done by plotting `log(abs(u))*sign(u)`"

# ╔═╡ 5b6d92a4-75d9-11eb-312d-25595fb8e8d5
let i = 2
	pl = PaynePolygon.plot_eigenfunction(
		domains[i], 
		us[i], 
		λs[i], 
		N,
		N,
		highlight_nodal_line = true, 
	)
	
	savefig(pl, "../figures/nodal-line.pdf")
	pl
end

# ╔═╡ Cell order:
# ╟─942342c0-705e-11eb-28c6-b17eb6e86323
# ╠═e800aea2-705e-11eb-37d1-c51771f6f70b
# ╠═11869a04-7062-11eb-1b79-dd45bea5b782
# ╠═f2626688-705e-11eb-02ab-0520ae1cbc31
# ╟─bbb33f6a-7769-11eb-1275-b74033697b16
# ╟─5a3367a6-7061-11eb-3c48-b937c8f9b006
# ╟─e1c8c5a2-776a-11eb-3a87-97f694a5429a
# ╟─767111ce-706a-11eb-3245-49883691e3d8
# ╟─1ace661c-7064-11eb-2b08-c92eafc00779
# ╟─764ff670-7066-11eb-259f-8b6a17f2f122
# ╟─0f902fe8-7068-11eb-3a97-914409e9fdf3
# ╟─8147efce-706b-11eb-3ca2-6b4edc35b3a4
# ╟─d354203a-706b-11eb-3378-5f6fe7f05596
# ╟─b01014d6-706a-11eb-1bd4-95c246347d8a
# ╟─f51966c8-776a-11eb-2667-f5671d29a85f
# ╟─5b6d92a4-75d9-11eb-312d-25595fb8e8d5
