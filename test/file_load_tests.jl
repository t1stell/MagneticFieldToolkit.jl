@testset "MagneticFieldToolkit - Load tests" begin
  mgridname = joinpath(@__DIR__, "mgrid_julia_test.nc")
  mg = read_mgrid(mgridname, [1.0E6])
  atol = 1.0E-8
  @testset "Load mgrid" begin
    @test mg.nfp == 4;
    @test size(mg.coords) == (10, 9, 15);
    @test size(mg.field_data[1]) == (10, 9, 15);
  end
  @testset "Mgrid interpolator" begin
    p1 = mg(1.41, 0.0, 0.0)
    @test isapprox(p1[1], 0.0, atol=atol)
    @test isapprox(p1[2], -86.910747255151, atol=atol)
    @test isapprox(p1[3], -49.432379093587, atol=atol)
    p2 = mg(1.42, 0.2, 0.1)
    @test isapprox(p2[1], 34.652657083547204, atol = atol)
    @test isapprox(p2[2], -77.66079279001298, atol = atol)
    @test isapprox(p2[3], -27.781203767316256, atol = atol)
  end
  @testset "Multiple current tests" begin
    mg = read_mgrid(mgridname, [1.0E6, -1.0E5, 0.0, 1.0E5])
    @test mg.nfp == 4;
    @test size(mg.coords) == (10, 9, 15);
    @test size(mg.field_data, 1) == 3;
    p1 = mg(1.41, 0.0, 0.0)
    @test isapprox(p1[1], 0.0, atol=atol)
    @test isapprox(p1[2], -86.7585154769, atol=atol)
    @test isapprox(p1[3], -49.3413864929, atol=atol)
    @test p1 == Tuple(compute_magnetic_field(mg, [1.41, 0.0, 0.0]))
    p2 = mg(1.42, 0.1, 0.1)
    @test isapprox(p2[1], 20.025683282642692, atol = atol)
    @test isapprox(p2[2], -78.886457262367, atol = atol)
    @test isapprox(p2[3], -41.46015776598917, atol = atol)
    @test p2 == Tuple(compute_magnetic_field(mg, [1.42, 0.1, 0.1]))
  end
  bmwname = joinpath(@__DIR__, "bmw_julia_test.nc")
  bmw = read_bmw(bmwname)
  @testset "Load bmw" begin
    @test bmw.nfp == 4;
    @test size(bmw.coords) == (10, 13, 11)
    @test size(bmw.field_data[1]) == (10, 13, 11)
    @test size(bmw.potential_data[1]) == (10, 13, 11)
  end
  p1 = bmw(2.0, 0.0, 0.0, A=true)
  p2 = bmw(2.1, 0.2, 0.1, A=true)
#  println("p1: ",p1)
#  println("p2: ",p2)
  
  @testset "BMW Interpolator" begin
    @test isapprox(p1[1][1], 0.0, atol=atol)
    @test isapprox(p1[1][2], 0.0013894885764634476, atol=atol)
    @test isapprox(p1[1][3], 0.002346906171789318, atol=atol)
    @test isapprox(p1[2][1], 0.0, atol=atol)
    @test isapprox(p1[2][2], 0.0008112768000118865, atol=atol)
    @test isapprox(p1[2][3], -5.2768618711445216e-5, atol=atol)
    @test isapprox(p2[1][1], 0.0020025045176464535, atol=atol)
    @test isapprox(p2[1][2], 0.002746327862181909, atol=atol)
    @test isapprox(p2[1][3], 0.003920916661087109, atol=atol)
    @test isapprox(p2[2][1], -0.00020813119536074314, atol=atol)
    @test isapprox(p2[2][2], 0.0014316070601660287, atol=atol)
    @test isapprox(p2[2][3], -6.159992591618763e-5, atol=atol)
  end
end
