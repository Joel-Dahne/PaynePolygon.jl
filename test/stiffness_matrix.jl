@testset "edges_to_triangles" begin
    edge_to_triangle1, edge_to_triangle2 = PaynePolygon.edges_to_triangles(36)
    # Compare to output from Matlab
    vars = matread("edges_to_triangles.mat")

    @test all(edge_to_triangle1 .== vars["EdgeToTriangle1"]')
    @test all(edge_to_triangle2 .== vars["EdgeToTriangle2"]')
end
