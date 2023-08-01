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
end


