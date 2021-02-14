@testset "bounded_by" begin
    RR = RealField(256)

    # TODO: Add some better tests for this
    @test PaynePolygon.bounded_by(sin, RR(0), 2RR(π), RR(1 + sqrt(eps())))
    @test PaynePolygon.bounded_by(sin, RR(0), 2RR(π), RR(1 + sqrt(eps())), use_taylor = true)
    @test !PaynePolygon.bounded_by(sin, RR(0), 2RR(π), RR(1 - sqrt(eps())))
    @test !PaynePolygon.bounded_by(sin, RR(0), 2RR(π), RR(1 - sqrt(eps())), use_taylor = true)
end
