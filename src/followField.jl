struct InterpolationParameters{T}
    itp::BFieldInterpolator{T}
    values::Vector{T}
    ϕ_max::T
    r_min::T
    r_max::T
    z_min::T
    z_max::T
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
  

function follow_field(itp::MagneticField{T},
                      rϕz::Array{Float64},
                      ϕ_end::Float64;
                      ϕ_step::Float64=zero(T),
                      poincare::Bool=false
                     ) where {T}
    ϕ_start = rϕz[2]
    u = @SVector [rϕz[1], rϕz[3]]
    ϕ_span = (ϕ_start,ϕ_end)
    params = InterpolationParameters(itp)
    prob = ODEProblem(field_deriv_ϕ, u, ϕ_span, params)
    if poincare
        N = abs(ϕ_end - ϕ_start)/(2π/itp.nfp)
        saveat = [i * 2*π/itp.nfp + ϕ_start for i in 1:N]
    else
        saveat = []
    end
    abs(ϕ_step) > zero(T) ? solve(prob, Tsit5(), dtmax = ϕ_step, saveat = saveat) :
                              solve(prob, Tsit5(), saveat = saveat)

end

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
end

function field_deriv_ϕ!(du::Vector{Float64},
                          u::Vector{Float64},
                          itp::BFieldInterpolator,
                          ϕ::Float64;
                         )
    r = u[1]
    z = u[2]
    br, bϕ, bz = itp(r, ϕ, z)
    du[1] = r * br/bϕ
    du[2] = r * bz/bϕ
end

function field_deriv_ϕ(u::AbstractVector,
                         itp::MagneticField{T},
                         ϕ::Float64;
                        ) where {T}
  ϕ = mod(ϕ, 2π/itp.nfp)
  br, bϕ, bz = itp(u[1], ϕ, u[2])
  dr = u[1] * br/bϕ
  dz = u[1] * bz/bϕ
  SVector{2,T}(dr, dz)
end

"""
"""
function field_deriv_ϕ(u::AbstractVector{T},
                         p::InterpolationParameters{T},
                         ϕ::T;
                        ) where {T}
    #ϕ = mod(ϕ, p.ϕ_max)
    if p.r_min < u[1] < p.r_max && p.z_min < u[2] < p.z_max
        p.values .= p.itp(u[1], ϕ, u[2])
        bϕ = u[1] / p.values[2]
        return SVector{2,T}(bϕ * p.values[1], bϕ * p.values[3])
    else
        return u
    end
end


#integration with respect to arclength
function field_deriv_s(u::AbstractVector,
                       p::InterpolationParameters{T},
                       s::Float64;) where {T}
    ϕ = mod(u[3], p.ϕ_max)
    map!(i->getfield(p.itp, i)(u[1], ϕ, u[2]), p.values, 1:3)
    bmagsq = sum(p.values.^2)
    dr = p.values[1]/bmagsq
    dz = p.values[2]/bmagsq
    dϕ = (p.values[3]/u[1])/bmagsq
    SVector{3, T}(dr, dϕ, dz)

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
