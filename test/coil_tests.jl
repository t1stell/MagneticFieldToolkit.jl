@testset "Coil tests" begin

  coilname = joinpath(@__DIR__, "circular_vmec_coils.txt")
  cset = read_vmec_coils(coilname)
  rtol = 1.0E-6
  @testset "verify coil reading" begin
    ncoils = 0
    @test length(cset.family) == 5
    for (i, family) in enumerate(cset.family)
      ncoils += length(family.coil)
      @test family.name == "coil"*string(i)
      for coil in family.coil
        @test length(coil.x) == 100
        @test length(coil.y) == 100
        @test length(coil.z) == 100
        @test coil.current == 1.0
      end
    end
  end
  @testset "verify basic calcs" begin
    @test isapprox(MagneticFieldToolkit.zextreme_coils(cset), 
                   0.9998741276738752, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.zextreme_coils(cset, zmax=false),
                   -0.9998741276738752, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.rextreme_coils(cset), 
                   11.000000000000002, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.rextreme_coils(cset, rmax=false), 
                   9.0005034576168122, rtol=rtol)
  end
end

