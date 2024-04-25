@testset "Field following tests" begin 
    #make a fake grid with field only in the toroidal direction
    rtol = 1.0E-10
    rtol_lo = 1.0E-5
    r = range(0.8, 1.2, 50)
    z = range(-0.2, 0.2, 50)
    θ = range(0, 2π/5, 100)
    fullSize = (length(r), length(θ), length(z))
    r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
    θ_grid = reshape(repeat(θ,inner=length(r),outer=length(z)),fullSize)
    z_grid = reshape(repeat(z,inner=length(r)*length(θ)),fullSize)
    field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))

    Br = zeros(length(r), length(θ), length(z))
    Bz = similar(Br)
    Bθ = similar(Br)
    Bz[:,:,:] .= 0.0
    Bθ[:,:,:] .= 1.0
    mg = MagneticField(field_coords, Br, Bθ, Bz, nfp=5)

    @testset "Follow a point in a perfectly toroidal field" begin
        for R in range(0.9, 1.1, 5)
            Z = (R-1) * 0.1 
            rθz = [R, 0.0, Z]
            a = follow_field(mg, rθz, 2π, ϕ_step = π/100, poincare = true, poincare_res = π/4)
            for k in a.u
                @test isapprox(k[1], R, rtol=rtol)
                @test isapprox(k[2], Z, rtol=rtol)
            end
        end
    end
    #make a wall that we can use for testing collisions
    nζ = 180
    nθ = 64
    θrange = range(0, 2π, nθ)
    ζs = range(0, π, nζ)
    Rp = Array{Float64}(undef, nθ, nζ)
    Zp = similar(Rp)
    θs = similar(Rp)
    for (iζ, ζ) in enumerate(ζs)
        a = 0.02+abs(π/2-ζ)/(π/2)*0.18
        Zp[:,iζ] = a .* sin.(θrange)
        Rp[:,iζ] = 1.0 .+ (a.* cos.(θrange))
        θs[:,iζ] = θrange
    end
    println(Rp[:,1])
    (R,Z) = StellaratorGrids.create_wall_splines(ζs, θrange, Rp, Zp, "")
    ves = FlareWall(nζ, nθ, Rp, Zp, collect(ζs), θs, 2, R, Z, "")
    m = π/2 / (1.02-1.2)
    b = -1.2*m
    ϕ_step = π/100
    @testset "Follow a point to a wall" begin
        for R in range(1.03, 1.18, 6)
            rθz = [R, 0.0, 0.0]
            (a, s) = follow_to_wall(mg, rθz, 2π, ves, false, ϕ_step = ϕ_step, rtol=1.0E-10)
            @test isapprox(R, a.u[end][1], rtol=rtol)
            #the toroidal direction is only valid to the resolution of the wall
            θg = m*R+b
            @test isapprox(a.t[end], θg, rtol=rtol_lo)
            @test s == 1
         end
    end
    @testset "Follow backwards" begin
        for R in range(1.03, 1.18, 6)
            rθz = [R, 0.0, 0.0]
            (a, s) = follow_to_wall(mg, rθz, -2π, ves, false, ϕ_step = ϕ_step, rtol=1.0E-10)
            @test isapprox(R, a.u[end][1], rtol=rtol)
            θg =  -(m*R+b)
            @test isapprox(a.t[end], θg, rtol=rtol_lo)
            @test s == 1
        end
    end

    @testset "Check out of bounds initial point" begin
        (a, s) = follow_to_wall(mg, [1.3, 0.2, 0.5], 2π, ves, false)
        @test isapprox(a.u[end][1], 1.3, rtol=rtol)
        @test isapprox(a.u[end][2], 0.5, rtol=rtol)
        @test isapprox(a.t[end], 0.2, rtol=rtol)
        @test s == 1
    end


    coilname = joinpath(@__DIR__, "circular_vmec_coils.txt")
    cset = read_vmec_coils(coilname)
    @testset "check following a coil set" begin
        a = follow_field(cset, [10.0, 0.0, 0.0], 2π, ϕ_step = π/100, poincare=true, poincare_res = π/4)
        for u in a.u
            @test isapprox(u[1], 10.0, rtol=rtol_lo)
            @test abs(u[2]) < 1.0E-13
        end
    end

    @testset "check diffuse function" begin
        cc = Cylindrical(1.03, 0.0, 0.0)
        ccn = MagneticFieldToolkit.diffuse(cc, mg, 0.1)
        xyz1 = CartesianFromCylindrical()(cc)
        xyz2 = CartesianFromCylindrical()(ccn)
        dist = sqrt(sum((xyz2 .- xyz1).^2))
        @test isapprox(dist, 0.1, rtol=rtol)
    end

    function make_fake_poincare(iota::Real, nfp::Integer)
        iota_per_period = iota/nfp
        R0 = 1.0
        a = 0.1
        npoints = 10000
        R_poincare = [R0 + a*cos(i*iota_per_period*2*π) for i in 1:npoints]
        Z_poincare = [a*sin(i*iota_per_period*2*π) for i in 1:npoints]
        return R_poincare, Z_poincare 
    end

    @testset "Test iota calculation" begin
        #test iota calculation with known iota
        nfp = 5
        target_iota = 1.2
        R_poincare, Z_poincare = make_fake_poincare(target_iota, nfp)
        iota = iota_from_poincare(R_poincare, Z_poincare, nfp)
        @test isapprox(target_iota, iota, rtol=rtol)#given resolutions, this should be exact
    end

    @testset "test_flux_calculation" begin
        
        target_iota = sqrt(2) #use an irrational value
        nfp = 5
        R_poincare, Z_poincare = make_fake_poincare(target_iota, nfp)
        
        itp = MagneticFieldToolkit.make_interpolation(R_poincare, Z_poincare)
        tflux = MagneticFieldToolkit.toroidal_flux_integration(mg, itp) 
        @test isapprox(tflux, π*0.1^2, rtol=rtol_lo)
    end



end


