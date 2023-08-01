@testset "Field following tests" begin 
    #make a fake grid with field only in the toroidal direction
    rtol = 1.0E-10
    r = range(0.8, 1.2, 50)
    z = range(-0.2, 0.2, 50)
    θ = range(0, π/5, 50)
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
    nζ = 80
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
    (R,Z) = StellaratorGrids.create_wall_splines(ζs, θrange, Rp, Zp, "")
    ves = FlareWall(nζ, nθ, Rp, Zp, collect(ζs), θs, 2, R, Z, "")
    m = π/2 / (1.02-1.2)
    b = -1.2*m
    ϕ_step = π/100
    @testset "Follow a point to a wall" begin
        for R in range(1.03, 1.18, 6)
            rθz = [R, 0.0, 0.0]
            a = follow_to_wall(mg, rθz, 2π, ves, false, ϕ_step = ϕ_step)
            @test isapprox(R, a.u[end][1], rtol=rtol)
            #the toroidal direction is only valid to the nearest pi/100
            θg = m*R+b
            @test a.t[end] > θg
            @test a.t[end] < θg + ϕ_step
         end
    end
    @testset "Follow backwards" begin
        for R in range(1.03, 1.18, 6)
            rθz = [R, 0.0, 0.0]
            a = follow_to_wall(mg, rθz, -2π, ves, false, ϕ_step = ϕ_step)
            @test isapprox(R, a.u[end][1], rtol=rtol)
            θg =  -(m*R+b)
            @test a.t[end] < θg
            @test a.t[end] > θg - ϕ_step 
        end
    end
end


