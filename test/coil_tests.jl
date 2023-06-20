@testset "Coil tests" begin

  coilname = joinpath(@__DIR__, "circular_vmec_coils.txt")
  cset = read_vmec_coils(coilname)
  rtol = 1.0E-6
  rtol_lo = 1.0E-4
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
    @test isapprox(MagneticFieldToolkit.extreme_coils(cset, :z), 
                   1.000545601674143, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.extreme_coils(cset, :z, vmax=false),
                   -1.0005456016741427, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.extreme_coils(cset, :r),
                   11.000000000000002, rtol=rtol)
    @test isapprox(MagneticFieldToolkit.extreme_coils(cset, :r, vmax=false), 
                   8.9998322371871672, rtol=rtol)
  end

  #create a set with just one coil for testing
  fil1 = cset.family[1].coil[1]
  fam1 = MagneticFieldToolkit.CoilFamily([fil1], "1coil");
  cs1coil = MagneticFieldToolkit.CoilSet([fam1])
  xyz = SVector(10.0, 0.0, 0.0)
  cc = Cylindrical(10.0, 0.0, 0.0)
  @testset "verify magnetic field on single current loop" begin
    @test isapprox(compute_magnetic_field(fil1, xyz), SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc), SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz), SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc), SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, xyz, current=10.0), SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc, current=10.0), SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz, currents=[10.0]), SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc, currents=[10.0]), SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
  end
  @testset "compute magnetic potential for single current loop" begin
    @test isapprox(compute_magnetic_potential(fil1, xyz), SVector(3.809452748892134e-22,  0.0,  2.0525360516352236e-12), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(fil1, cc), SVector(3.809452748892134e-22,  0.0,  2.0525360516352236e-12), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cs1coil, xyz), SVector(3.809452748892134e-22,  0.0,  2.0525360516352236e-12), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cs1coil, cc), SVector(3.809452748892134e-22,  0.0,  2.0525360516352236e-12), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(fil1, xyz, current=10.0), SVector(3.809452748892134e-21,  0.0,  2.0525360516352236e-11), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(fil1, cc, current=10.0), SVector(3.809452748892134e-21,  0.0,  2.0525360516352236e-11), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cs1coil, xyz, currents=[10.0]), SVector(3.809452748892134e-21,  0.0,  2.0525360516352236e-11), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cs1coil, cc, currents=[10.0]), SVector(3.809452748892134e-21,  0.0,  2.0525360516352236e-11), rtol = rtol_lo)
  end

  #test the individual coil loader
  @testset "load individual coils" begin
    fs1 = joinpath(@__DIR__, "sample_coil_1.txt")
    fs2 = joinpath(@__DIR__, "sample_coil_2.txt")
    cs = read_coil_files([fs1, fs2], 5)
    @test length(cs) == 10
    @test isapprox(cs[1].x(0), 11.0, rtol=rtol)
    @test isapprox(cs[2].x(0), 10.947031993394164, rtol=rtol)
    @test isapprox(cs[3].x(0), 11.0 * cos(2π/5), rtol=rtol)
    @test isapprox(cs[5].x(0), 11.0 * cos(4π/5), rtol=rtol)
    @test isapprox(cs[7].x(0), 11.0 * cos(6π/5), rtol=rtol)
    @test isapprox(cs[9].x(0), 11.0 * cos(8π/5), rtol=rtol)
  end
end

