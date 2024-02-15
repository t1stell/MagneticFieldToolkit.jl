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
    @test isapprox(compute_magnetic_field(fil1, xyz, exact=true), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc, exact=true), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz, exact=true),
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc, exact=true), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, xyz, current=10.0, exact=true), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc, current=10.0, exact=true), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz, currents=[10.0], exact=true), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc, currents=[10.0], exact=true), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
  end
  @testset "do simplified biot-savart calculation" begin
    @test isapprox(compute_magnetic_field(fil1, xyz), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz),
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc), 
          SVector(0.0, -2*π*1.0E-7, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, xyz, current=10.0), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(fil1, cc, current=10.0), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, xyz, currents=[10.0]), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
    @test isapprox(compute_magnetic_field(cs1coil, cc, currents=[10.0]), 
          SVector(0.0, -2*π*1.0E-6, 0.0), rtol = rtol_lo)
  end

  #test the individual coil loader
  @testset "load individual coils" begin
    fs1 = joinpath(@__DIR__, "sample_coil_1.txt")
    fs2 = joinpath(@__DIR__, "sample_coil_2.txt")
    cs = read_coil_files([fs1, fs2], 5)
    @test length(cs.family[1].coil) == 10
    @test isapprox(cs.family[1].coil[1].x(0), 11.0, rtol=rtol)
    @test isapprox(cs.family[1].coil[2].x(0), 10.947031993394164, rtol=rtol)
    @test isapprox(cs.family[1].coil[3].x(0), 11.0 * cos(2π/5), rtol=rtol)
    @test isapprox(cs.family[1].coil[5].x(0), 11.0 * cos(4π/5), rtol=rtol)
    @test isapprox(cs.family[1].coil[7].x(0), 11.0 * cos(6π/5), rtol=rtol)
    @test isapprox(cs.family[1].coil[9].x(0), 11.0 * cos(8π/5), rtol=rtol)
  end
  
  @testset "compute magnetic potential for coil set, tested against mgrid" begin
    currents = [10.0, 10.0, 10.0, 10.0, 10.0]
    @test isapprox(compute_magnetic_potential(cset, xyz), 
          SVector(0.0,  0.0,  -1.0816251E-7), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, cc), 
          SVector(0.0,  0.0,  -1.0816251E-7), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, xyz, currents=currents), 
          SVector(0.0,  0.0,  -1.0816251E-6), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, cc, currents=currents), 
          SVector(0.0,  0.0,  -1.0816251E-6), rtol = rtol_lo)
  end
  @testset "compute magnetic potential for coil set (exact), tested against mgrid" begin
    currents = [10.0, 10.0, 10.0, 10.0, 10.0]
    @test isapprox(compute_magnetic_potential(cset, xyz, exact=true), 
          SVector(0.0,  0.0,  -1.0815561E-7), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, cc, exact=true), 
          SVector(0.0,  0.0,  -1.0815561E-7), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, xyz, currents=currents, exact=true), 
          SVector(0.0,  0.0,  -1.0815561E-6), rtol = rtol_lo)
    @test isapprox(compute_magnetic_potential(cset, cc, currents=currents, exact=true), 
          SVector(0.0,  0.0,  -1.0815561E-6), rtol = rtol_lo)
  end

  #create an mgrid from the single coil set, note this is note enough resolution
  #for an accurate answer.  more
  mg1coil = generate_mgrid(cs1coil, 20, 30, 20, 1)
  (bt, at) = mg1coil(10.2, 0.2, 0.2, [1.0E7])
  @testset "create mgrid from single coil" begin
    (b, a) = mg1coil(10.0, 0.0, 0.0, [1.0])
    @test abs(b[1]) < 1.0E-20  
    @test abs(b[3]) < 1.0E-20  
    @test abs(b[2]) < 6.4E-7 && abs(b[2]) > 6.2E-7 #converges to 2πE-7 at higher res
  end
  #test saving and reading a file
  savename = joinpath(@__DIR__, "single_coil_mgrid.h5")
  MagneticFieldToolkit.save_mgrid(mg1coil, savename)  
  @testset "writing/reading an mgrid h5 file" begin
    mgh5 = read_mgrid_h5(savename)
    (b, a) = mgh5(10.0, 0.0, 0.0, [1.0])

    @test abs(b[1]) < 1.0E-20  
    @test abs(b[3]) < 1.0E-20  
    @test abs(b[2]) < 6.4E-7 && abs(b[2]) > 6.2E-7 #converges to 2πE-7 at higher res
    (bt2, at2) = mgh5(10.2, 0.2, 0.2, [1.0e7])
    @test isapprox(bt[1], bt2[1], rtol=rtol_lo)
    @test isapprox(bt[2], bt2[2], rtol=rtol_lo)
    @test isapprox(bt[3], bt2[3], rtol=rtol_lo)
    @test isapprox(at[1], at2[1], rtol=rtol_lo)
    @test isapprox(at[2], at2[2], rtol=rtol_lo)
    @test isapprox(at[3], at2[3], rtol=rtol_lo)
    
    #test the other way of loading
    mgh5 = read_mgrid_h5(savename, currents=[1.0e7])
    (bt2, at2) = mgh5(10.2, 0.2, 0.2, A=true)
    @test isapprox(bt[1], bt2[1], rtol=rtol_lo)
    @test isapprox(bt[2], bt2[2], rtol=rtol_lo)
    @test isapprox(bt[3], bt2[3], rtol=rtol_lo)
    @test isapprox(at[1], at2[1], rtol=rtol_lo)
    @test isapprox(at[2], at2[2], rtol=rtol_lo)
    @test isapprox(at[3], at2[3], rtol=rtol_lo)


  end 
  rm(savename)

end

