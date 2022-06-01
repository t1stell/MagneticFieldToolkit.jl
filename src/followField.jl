struct InterpolationParameters{T}
    itp::BFieldInterpolator{T}
    values::Vector{T}
    ϕ_max::T
end

function follow_field(itp::BFieldInterpolator{T},
                      rzp::Array{Float64},
                      phi_end::Float64;
                      phi_step::Float64=zero(T),
                      poincare::Bool=false
                     ) where {T}
    phiStart = rzp[3]
    u = @SVector [rzp[1], rzp[2]]
    phi_span = (phiStart,phi_end)
    params = InterpolationParameters(itp, zeros(T,3), 2π/itp.nfp)
    prob = ODEProblem(field_deriv_phi, u, phi_span, params)
    if poincare
        N = abs(phi_end - phiStart)/(2π/itp.nfp)
        saveat = [i * 2*π/itp.nfp + phiStart for i in 1:N]
    else
        saveat = []
    end
    abs(phi_step) > zero(T) ? solve(prob, Tsit5(), dtmax = phi_step, saveat = saveat) :
                              solve(prob, Tsit5(), saveat = saveat)

end

function follow_field_s(itp::BFieldInterpolator{T},
                      rzp::Array{Float64},
                      s_end::Float64;
                      s_step::Float64=zero(T),
                      ) where{T}
    u = @SVector [rzp[1], rzp[2], rzp[3]]
    sSpan = (0, s_end)
    params = InterpolationParameters(itp, zeros(T,3), 2π/itp.nfp)
    prob = ODEProblem(field_deriv_s, u, sSpan, params)
    abs(s_step) > zero(T) ? solve(prob, Tsit5(), dtmax = s_step, saveat = s_step) :
                            solve(prob, Tsit5(), saveat = s_step)
end

function field_deriv_phi!(du::Vector{Float64},
                          u::Vector{Float64},
                          itp::BFieldInterpolator,
                          ϕ::Float64;
                         )
    r = u[1]
    z = u[2]
    br, bz, bp = itp(r, z, ϕ)
    du[1] = r * br/bp
    du[2] = r * bz/bp
end

function field_deriv_phi(u::AbstractVector,
                         itp::BFieldInterpolator{T},
                         ϕ::Float64;
                        ) where {T}
  ϕ = mod(ϕ, 2π/itp.nfp)
  br, bz, bp = itp(u[1], u[2], ϕ)
  dr = u[1] * br/bp
  dz = u[1] * bz/bp
  SVector{2,T}(dr, dz)
end

"""
"""
function field_deriv_phi(u::AbstractVector{T},
                         p::InterpolationParameters{T},
                         ϕ::T;
                        ) where {T}
    ϕ = mod(ϕ, p.ϕ_max)
    map!(i->getfield(p.itp, i)(u[1], u[2], ϕ), p.values, 1:3)
    bp = u[1] / p.values[3]
    SVector{2,T}(bp * p.values[1], bp * p.values[2])
end


#integration with respect to arclength
function field_deriv_s(u::AbstractVector,
                       p::InterpolationParameters{T},
                       s::Float64;) where {T}
    ϕ = mod(u[3], p.ϕ_max)
    map!(i->getfield(p.itp, i)(u[1], u[2], ϕ), p.values, 1:3)
    bmagsq = sum(p.values.^2)
    dr = p.values[1]/bmagsq
    dz = p.values[2]/bmagsq
    dϕ = (p.values[3]/u[1])/bmagsq
    SVector{3, T}(dr, dz, dϕ)

end

function poincare(itp::BFieldInterpolator{T},
                  initial_conditions::AbstractArray;
                  trace_ntransits::Integer = 1,
                  trace_nfp::Integer = 0,
                  ϕ_step::T=zero(T),
                 ) where {T}
  ϕ_nfp = 2π/itp.nfp
  N = iszero(trace_nfp) ? 2π * trace_ntransits : ϕ_nfp * trace_nfp
  ϕ_start = last(first(initial_conditions))
  ϕ_end = N * ϕ_nfp + ϕ_start 
  ϕ_span = (ϕ_start, ϕ_end)
  params = InterpolationParameters(itp, zeros(T,3), ϕ_nfp)
  u0 = SVector{2,T}(first(initial_conditions)[1:2])
  ϕ_saveat = [i * ϕ_nfp + ϕ_start for i in 1:N]
  prob = iszero(ϕ_step) ? ODEProblem(field_deriv_phi, u0, ϕ_span, params, saveat = ϕ_saveat) :
                          ODEProblem(field_deriv_phi, u0, ϕ_span, params, saveat = ϕ_saveat, dtmax = ϕ_step)

  function prob_func(prob, i, repeat)
    @info "Remaking at iteration $i with initial condition $(initial_conditions[i][1:2])"
    remake(prob, u0 = SVector{2,T}(initial_conditions[i][1:2]))
  end

  poincare_prob = EnsembleProblem(prob, prob_func = prob_func)
  solve(poincare_prob, Tsit5(), EnsembleThreads(), trajectories = length(initial_conditions))
end

function poincare(bfield::BField,
                  r₀::Union{T, AbstractVector{T}},
                  z₀::Union{T, AbstractVector{T}},
                  ϕₒ::T;
                  trace_ntransits::Integer = 1,
                  trace_nfp::Integer = 0,
                  ϕ_step::T=zero(T),
                 ) where {T}
    itp = BFieldInterpolator(bfield)
    full_size = (length(r₀),length(z₀))
    r_grid = reshape(repeat(rₒ, inner= length(zₒ)), full_size)
    z_grid = reshape(repeat(zₒ, outer= length(rₒ)), full_size)
    init_cond = Matrix{T}(undef, full_size)
    for i in eachindex(r_grid, z_grid, init_cond)
        init_cond[i] = [r_grid[i], z_grid[i], ϕ₀]
    end
    poincare(itp, init_cond; trace_ntransits = trace_ntransits, trace_nfp = trace_nfp, ϕ_step = ϕ_step)
end
