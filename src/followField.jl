struct InterpolationParameters{T}
    field_info::Union{MagneticField{T}, CoilSet{T}}
    values::Vector{T}
    ϕ_max::T
    r_min::T
    r_max::T
    z_min::T
    z_max::T
end

struct CoilInterpolationParameters{T}
    cs::CoilSet{T}
    values::Vector{T}
end

"""
  Constructor for InterpolationParameters Struct
"""
function InterpolationParameters(itp::MagneticField{T}) where{T}
  values = zeros(T, 3)
  ϕ_max = maximum(itp.coords.θ)
  r_max = maximum(itp.coords.r)
  r_min = minimum(itp.coords.r)
  z_max = maximum(itp.coords.z)
  z_min = minimum(itp.coords.z)
  return InterpolationParameters{T}(itp, values, ϕ_max, r_min, r_max, z_min, z_max)
end

function InterpolationParameters(cset::CoilSet{T}) where {T}
  values = zeros(T, 3)
  ϕ_max = 2π
  r_max = extreme_coils(cset, :r, vmax=true)
  r_min = extreme_coils(cset, :r, vmax=false)
  z_max = extreme_coils(cset, :z, vmax=true)
  z_min = extreme_coils(cset, :z, vmax=false)
  return InterpolationParameters{T}(cset, values, ϕ_max, r_min, r_max, z_min, z_max)
end

function diffuse(cc::Cylindrical, fieldinfo::Union{MagneticField{T}, CoilSet{T}}, D::Float64) where {T}
    xyz = CartesianFromCylindrical()(cc)
    rand_vec = randn(3) #random isotropic vector
    Bxyz = compute_magnetic_field(fieldinfo, xyz)
    rand_perp = cross(rand_vec, Bxyz) #random perpendicular vector
    rand_perp /= norm(rand_perp)
    rand_perp *= D
    xyzn = xyz .+ rand_perp
    return CylindricalFromCartesian()(xyzn)
end
      
  
function follow_field(fieldinfo::Union{MagneticField{T}, CoilSet{T}},
                      rϕz::Array{Float64},
                      ϕ_end::Float64;
                      ϕ_step::Float64=π/25,
                      poincare::Bool=false,
                      poincare_res::Real=2π,
                      maxiters::Int = 10^7,
                      include_first::Bool = true,
                      include_last::Bool = true
                     ) where {T}
    ϕ_start = rϕz[2]
    u = @MVector [rϕz[1], rϕz[3]]

    if include_first
        start_index = 0
    else
        start_index = 1
    end

    if include_last
        final_adjust = 0
    else
        final_adjust = 1
    end

    #if poincare
    #    ϕ_end = ϕ_end + poincare_res / 2.0
    #end
    ϕ_span = (ϕ_start,ϕ_end)
    params = InterpolationParameters(fieldinfo)
    prob = ODEProblem(field_deriv_ϕ, u, ϕ_span, params,maxiters=maxiters)
    if poincare
        if poincare_res == 2π
          ϕ_max = params.ϕ_max
        else
          ϕ_max = poincare_res
        end
        N = abs(ϕ_end - ϕ_start)/(ϕ_max)
        saveat = [i * (ϕ_max) + ϕ_start for i in start_index:N - final_adjust]
    else
        saveat = []
    end
    abs(ϕ_step) > zero(T) ? solve(prob, Tsit5(), dtmax = ϕ_step, saveat = saveat) :
                              solve(prob, Tsit5(), saveat = saveat)

end
"""
    function follow_to_wall(fieldinfo, rϕz, ϕ_end, wall, wall_inverse;
                            ϕ_step = π/25, poincare=false,
                            poincare_res = 2π, wall_res = 128, rtol = 1.0E-6
                            diffusion = 0.0)

This function is like follow_field except it carries around a bunch of extra stuff so it can
also calculate intersections with the wall. 

#Arguments

 - `field_info::Union{MagneticField{T}, CoilSet{T}}`: The magnetic field, it can either be a magnetic field object (i.e. Br, Bz, Bϕ on a cylindrical grid) or a CoilSet object. If speed is necessary, should try to generate a magnetic field object
 - `rϕz::Array{Float64}`: The starting point in r, ϕ, z where ϕ is the cylindrical coordinate. (be careful, θ is used as the cylindrical coordinate in CoordinateTransformations)
 - `ϕ_end::Float64`: The value at which to terminate the following.  If less than the ϕ value in `rϕz`, it will follow in the reverse direction
 - `wall::FlareWall`: Can be either a singular wall, or an array of walls.  These objects are found in StellaratorGrids.jl
 - `wall_inverse::Bool`: Can be a single value or an array of the same size as FlareWall. If set to "false" points outisde the wall polygon will be marked as out of bounds.  If true, points inside the wall will be marked out of bounds. 

#Optional Arguments
 - `ϕ_step`: default π/25.  This is the internal integration steps. Higher values will calculate quicker, but may produce unsatisfactory errors.  Each configuration will need to do a convergence scan to determine the best value
 - `poincare`: This is a bool value to output a poincare result (i.e. export every value at a given phi)
 - `poincare_res`: default 2π.  This is the value at which to output a poincare result.  
 - `wall_res`: default 128.  The value on which to resample the wall for the purposes of quick polygonal calculations
 - `rtol`: default 1.0E-6. This is the tolerance for the line search segment to find the wall intersection. At low values, the proper resolution is actually set by the wall resolution
 - `diffusion`: default 0.  This is the value in meters^2/meters. That is, every meter step perpendicularly this value. The actual diffusion is calculated at every integration step
"""
function follow_to_wall(fieldinfo::Union{MagneticField{T}, CoilSet{T}},
                        rϕz::Array{Float64},
                        ϕ_end::Real,
                        wall::Vector{FlareWall},
                        wall_inverse::Vector{Bool};
                        ϕ_step::Real=π/25,
                        poincare::Bool=false,
                        poincare_res::Real=2π,
                        wall_res::Integer = 128, #eventually allow this to be a vector
                        rtol::Float64=1.0E-6, #tolerance for the last step
                        diffusion::Float64=0.0,
                        maxiters::Int = 10^7
                       ) where {T}
    rwall = Vector{T}(undef, wall_res)
    zwall = Vector{T}(undef, wall_res)
    ϕ_start = rϕz[2] #starting toroidal angle
    u = @MVector [rϕz[1], rϕz[3]] #starting u vector, (R, Z).  
    start_inside = true #flag for whether we start inside the boundary
    last_good = rϕz[:] #tracker for the last good point
    penult_good = rϕz[:] #tracker for the penultimate point (needed for diffusion)
    direction = sign(ϕ_end - ϕ_start)

    # this is a function used in the callback, it checks if we are currently on the "correct" side
    # of all the walls
    function outside_bounds(u::MVector{2, F}, t::F) where {F <: AbstractFloat}
        #println("checking bounds") 
        outside = 0
        θs = range(0, 2π, wall_res)

        #different work flows depending on the type of wall (TODO: move these subfunctions into StellaratorGrids and use multiple dispatch)
        for i in 1:length(wall)
            ζmax = 2*π/wall[i].nfp
            ζ = mod(t, ζmax)
            #note the wall is periodic by default, will need to fix that for finite extent walls
            if ζ < wall[i].ζs[1] || ζ > wall[i].ζs[end]
                continue
            end
            for (j,θ) in enumerate(θs)
                rwall[j] = wall[i].R(θ, ζ)
                zwall[j] = wall[i].Z(θ, ζ)
            end
            if !in_surface(u[1], u[2], rwall, zwall, inverse=wall_inverse[i])
                return i
            end
        end
        penult_good = copy(last_good) #if we don't copy, it will just pass by reference
        last_good = [u[1], t, u[2]]
        return outside
    end

    function outside_bounds_bool(u::MVector{2, F}, t::F) where {F <: AbstractFloat}
        v = outside_bounds(u, t)
        if v == 0
            return false
        else
            return true
        end
    end
    
    function diffuse_affect!(integrator::OrdinaryDiffEq.ODEIntegrator)
        #println("diffusing")
        if diffusion <= 0.0
            return
        end
        u = integrator.u
        t = integrator.t
        cc = Cylindrical(u[1], t, u[2])
        dist = sqrt(sum((last_good .- penult_good).^2))
        ccn = diffuse(cc, fieldinfo, dist*diffusion)
        #offset = round((cc.θ - ccn.θ)/π) * π
        #println(cc.θ," ",ccn.θ," ",offset)

        integrator.u[1] = ccn.r
        integrator.u[2] = ccn.z

        #Updating the θ value seemed to cause a memory leak for unknown reasons
        #we'll comment this out for now
        #integrator.t = ccn.θ + offset

        if (direction < 0 && integrator.t < ϕ_end) || (direction > 0 && integrator.t > ϕ_end)
            #integrator.t = ϕ_end
            terminate!(integrator)
        end
       
        #this check slows things down considerably, let's remove it 
        #if outside_bounds_bool(integrator.u, integrator.t)
        #    terminate!(integrator)
        #end
        
    end
    
    #before starting check if we're outside bounds
    #don't return immediately.  To retain the correct output structure, make the follower immediately exit
    #by setting end = start
    if outside_bounds_bool(u, ϕ_start)
        ϕ_end = ϕ_start
        struck_target = outside_bounds(u, ϕ_start)
        start_inside = false
    end


    ϕ_span = (ϕ_start,ϕ_end)
    params = InterpolationParameters(fieldinfo) #struct with calculation info
    prob = ODEProblem(field_deriv_ϕ, u, ϕ_span, params, maxiters=maxiters) #the ODE problem to solve
    condition_wall(u, t, integrator) = outside_bounds_bool(u, t) #The bound condition
    stop_affect!(integrator) = terminate!(integrator) #identify the affect used for the bound condition
    condition_always(u, t, integrator) = true
    cb_wall = DiscreteCallback(condition_wall, stop_affect!) #set up the actual bound
    #diffusion goes here
    cb_diffusion = DiscreteCallback(condition_always, diffuse_affect!) # do this at every integration step
    cbset = CallbackSet(cb_diffusion, cb_wall) #make it into a set, both are DiscreteCallbacks so they'll be checked in order

    #This section creates an array of values for the integrator to save at, at the resolution of poincare_res
    if poincare
        if poincare_res == 2π
          ϕ_max = params.ϕ_max
        else
          ϕ_max = poincare_res
        end
        N = abs(ϕ_end - ϕ_start)/(ϕ_max)
        saveat = [i * (ϕ_max) + ϕ_start for i in 1:N]
    else
        saveat = []
    end
    a = abs(ϕ_step) > zero(T) ? solve(prob, Tsit5(), dtmax = ϕ_step, saveat = saveat, callback=cbset) :
                              solve(prob, Tsit5(), saveat = saveat, callback=cbset)

    #Find the last step, we do this by starting at the last good point and gradually shrinking
    if !start_inside
        return a, struck_target #we checked this earlier
    end

    tolcheck = abs(a.t[end] - last_good[2]) #distance between the integrator final point, and the last known point inside the boundary
    dϕ_new = tolcheck/10 #New integrator step
    c = 0
    t_end = a.t[end]
    u_end = a.u[end]
    #find the target
    struck_target = outside_bounds(u_end, t_end) 
    while tolcheck > rtol
        u = @MVector [last_good[1], last_good[3]] 
        ϕ_span = (last_good[2], a.t[end])
        prob = ODEProblem(field_deriv_ϕ, u, ϕ_span, params, maxiters=maxiters)
        #note do not use diffusion for this, do not use the full cbset
        a_new = solve(prob, Tsit5(), dtmax = dϕ_new, callback = cb_wall) #set up a new problem
        t_end = a_new.t[end]
        u_end = a_new.u[end]
        tolcheck = abs(a_new.t[end] - last_good[2]) #new difference between end and last good point
        dϕ_new = tolcheck/10 #divide integrator step by another factor of 10 and repeat
        c += 1
        if c == 15 # This is a check, at this point, you are probably at machine precision.  Should never get here
            break
        end
    end

    # add the points to the end (note there will be some bad points prior to this, should we remove them?`
    push!(a.t, t_end) 
    push!(a.u, u_end)
    return a, struck_target
end

function follow_to_wall(fieldinfo::Union{MagneticField{T}, CoilSet{T}},
                        rϕz::Array{Float64},
                        ϕ_end::Real,
                        wall::FlareWall,
                        wall_inverse::Bool;
                        ϕ_step::Real=π/25,
                        poincare::Bool=false,
                        poincare_res::Real=2π,
                        wall_res::Integer = 128, #eventually allow this to be a vector
                        rtol::Float64=1.0E-6, #tolerance for the last step
                        diffusion::Float64=0.0,
                        maxiters::Int = 10^7
                        ) where {T}

    wall_list = Vector{FlareWall}(undef, 1)
    wall_list[1] = wall
    return follow_to_wall(fieldinfo, rϕz, ϕ_end, wall_list, [wall_inverse], 
                          ϕ_step = ϕ_step, poincare=poincare, poincare_res=poincare_res,
                          wall_res = wall_res, rtol=rtol, diffusion=diffusion, maxiters=maxiters)
end

"""
function follow_field_s(itp::MagneticField{T},
                      rϕz::Array{Float64},
                      s_end::Float64;
                      s_step::Float64=zero(T),
                      ) where{T}
    u = @SVector [rϕz[1], rϕz[2], rϕz[3]]
    s_span = (0, s_end)
    params = InterpolationParameters(itp, zeros(T,3), 2π/itp.nfp)
    prob = ODEProblem(field_deriv_s, u, s_span, params)
    abs(s_step) > zero(T) ? solve(prob, Tsit5(), dtmax = s_step, saveat = s_step) :
                            solve(prob, Tsit5(), saveat = s_step)
"""
#This version gets called from the main follow field
#it picks out the interpolator type via the field_info field
#and passes to the correct field_deriv_ϕ function below
function field_deriv_ϕ(u::AbstractVector{T},
                         p::InterpolationParameters{T},
                         ϕ::T;
                        ) where {T}
    if p.r_min < u[1] < p.r_max && p.z_min < u[2] < p.z_max
        return field_deriv_ϕ(u, p.field_info, ϕ)
    else
        return u
    end
end


function field_deriv_ϕ(u::AbstractVector,
                         itp::MagneticField{T},
                         ϕ::Float64;
                        ) where {T}
  ϕ = mod(ϕ, 2π/itp.nfp)
  cc = Cylindrical(u[1], ϕ, u[2])
  brϕz = compute_magnetic_field(itp, cc) 
  dr = u[1] * brϕz[1]/brϕz[2]
  dz = u[1] * brϕz[3]/brϕz[2]
  return MVector{2,T}(dr, dz)
end

function field_deriv_ϕ( u::AbstractVector,
                         cset::CoilSet{T},
                         ϕ::Float64;
                        ) where {T}
    
    r = u[1]
    z = u[2]
#    cc = Cylindrical(r, ϕ, z)
    (br, bϕ, bz) = compute_magnetic_field(cset, [r, ϕ, z])
    dr = r * br/bϕ
    dz = r * bz/bϕ
    return MVector{2,T}(dr, dz)
end

#todo: fix this. right now it doesn't seem to know how far to follow
function follow_field_s(fieldinfo::Union{MagneticField{T}, CoilSet{T}},
                      rϕz::Array{Float64},
                      s_end::Float64;
                      s_step::Float64=zero(T),
                      ) where{T}
    u = @SVector [rϕz[1], rϕz[2], rϕz[3]]
    s_span = (0, s_end)
    params = InterpolationParameters(fieldinfo)
    prob = ODEProblem(field_deriv_s, u, s_span, params)
    saveat = []
    abs(s_step) > zero(T) ? solve(prob, Tsit5(), dtmax = s_step, saveat = saveat) :
                            solve(prob, Tsit5(), saveat = saveat)
end

#Check bounds and then call the correct derivative function
function field_deriv_s(u::AbstractVector{T},
                         p::InterpolationParameters{T},
                         s::T;
                        ) where {T}
    if p.r_min < u[1] < p.r_max && p.z_min < u[2] < p.z_max
        return field_deriv_s(u, p.field_info, s)
    else
        return u
    end
end

function field_deriv_ϕ(u::AbstractVector{T},
                         p::CoilInterpolationParameters{T},
                         ϕ::T;
                        ) where {T}
     #ϕ = mod(ϕ, p.ϕ_max)
     p.values .= compute_magnetic_field(p.cs, Cylindrical(u[1], ϕ, u[2]))
     #println(p.values, u, ϕ)
     bϕ = u[1] / p.values[2]
     return SVector{2,T}(bϕ * p.values[1], bϕ * p.values[3])
end

#integration with respect to arclength
function field_deriv_s(u::AbstractVector,
                       itp::MagneticField{T},
                       s::Float64;) where {T}
    ϕ = mod(u[2], 2π/itp.nfp)
    br, bϕ, bz = itp(u[1], ϕ, u[3])
    bmagsq = br^2 + bϕ^2 + bz^2
    dr = br/bmagsq
    dz = bz/bmagsq
    dϕ = (bϕ/u[1])/bmagsq
    return SVector{3, T}(dr, dϕ, dz)
end

function poincare(itp::MagneticField,
                  initial_conditions::AbstractArray{Vector{T}};
                  trace_ntransits::Integer = 1,
                  trace_nfp::Integer = 0,
                  ϕ_step::T=zero(T),
                  maxiters::Int = 10^5,
                 ) where {T}
  ϕ_nfp = 2π/itp.nfp
  N = iszero(trace_nfp) ? 2π * trace_ntransits : ϕ_nfp * trace_nfp
  ϕ_start = first(initial_conditions)[2]
  ϕ_end = N * ϕ_nfp + ϕ_start 
  ϕ_span = (ϕ_start, ϕ_end)
  params = InterpolationParameters(itp)
  u0 = SVector{2,T}(first(initial_conditions)[1],first(initial_conditions)[3])
  ϕ_saveat = [i * ϕ_nfp + ϕ_start for i in 1:N]
  prob = iszero(ϕ_step) ? ODEProblem(field_deriv_ϕ, u0, ϕ_span, params, saveat = ϕ_saveat, maxiters = maxiters) :
                          ODEProblem(field_deriv_ϕ, u0, ϕ_span, params, saveat = ϕ_saveat, maxiters = maxiters, dtmax = ϕ_step)

  function prob_func(prob, i, repeat)
    @debug "Remaking trajectory $i with initial condition ($(initial_conditions[i][1]), $(initial_conditions[i][3])))"
    remake(prob, u0 = SVector{2,T}(initial_conditions[i][1], initial_conditions[i][3]))
  end

  poincare_prob = EnsembleProblem(prob, prob_func = prob_func)
  solve(poincare_prob,
        Tsit5(),
        EnsembleThreads(),
        trajectories = length(initial_conditions), 
       )
end

function poincare(itp::MagneticField,
                  r₀::Union{T, AbstractVector{T}},
                  z₀::Union{T, AbstractVector{T}},
                  ϕ₀::T;
                  trace_ntransits::Integer = 1,
                  trace_nfp::Integer = 0,
                  ϕ_step::T=zero(T),
                  maxiters::Int = 10^5,
                 ) where {T}
    itp = BFieldInterpolator(bfield)
    full_size = (length(r₀),length(z₀))
    r_grid = reshape(repeat(r₀, inner= length(z₀)), full_size)
    z_grid = reshape(repeat(z₀, outer= length(r₀)), full_size)
    init_cond = Matrix{Vector{T}}(undef, full_size)
    for i in eachindex(r_grid, z_grid, init_cond)
        init_cond[i] = [r_grid[i], ϕ₀, z_grid[i]]
    end

    poincare(itp, init_cond;
             trace_ntransits = trace_ntransits,
             trace_nfp = trace_nfp,
             ϕ_step = ϕ_step,
             maxiters = maxiters)
end
