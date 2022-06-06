@testset "MagneticFieldToolkit - Load tests" begin
  mgridname = joinpath(@__DIR__, "mgrid_julia_test.nc")
  mg = read_mgrid(mgridname, [1.0E6])
  atol = 1.0E-5
  @testset "Load mgrid" begin
    @test mg.nfp == 4;
    @test size(mg.coords) == (10, 15, 8);
    @test size(mg.field_data, 1) == 3;
  end
  @testset "Mgrid interpolator" begin
    p1 = mg(1.41, 0.0, 0.0)
    @test isapprox(p1[1], 0.0, atol=atol)
    @test isapprox(p1[2], -49.432379093587, atol=atol)
    @test isapprox(p1[3], -86.910747255151, atol=atol)
    p2 = mg(1.42, 0.1, 0.1)
    @test isapprox(p2[1], 24.83161930676, atol = atol)
    @test isapprox(p2[2], -41.6860661569, atol = atol)
    @test isapprox(p2[3], -79.6157100469, atol = atol)
  end
  @testset "Multiple current tests" begin
    mg = read_mgrid(mgridname, [1.0E6, -1.0E5, 0.0, 1.0E5])
    @test mg.nfp == 4;
    @test size(mg.coords) == (10, 15, 8);
    @test size(mg.field_data, 1) == 3;
    p1 = mg(1.41, 0.0, 0.0)
    @test isapprox(p1[1], 0.0, atol=atol)
    @test isapprox(p1[2], -49.3413864929, atol=atol)
    @test isapprox(p1[3], -86.7585154769, atol=atol)
    p2 = mg(1.42, 0.1, 0.1)
    @test isapprox(p2[1], 24.85618512831, atol = atol)
    @test isapprox(p2[2], -41.6001501965, atol = atol)
    @test isapprox(p2[3], -79.4934865911, atol = atol)


  end
end
