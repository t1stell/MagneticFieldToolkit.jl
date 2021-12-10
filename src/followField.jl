struct InterpolationParameters{T}
  itp::BFieldInterpolator{T}
  values::Vector{T}
  ϕ_max::T
end

function followField(itp::BFieldInterpolator{T},
                     rzp::Array{Float64},
                     phiEnd::Float64;
                     phiStep::Float64=zero(T),
                     poincare::Bool=false
                    ) where {T}
  phiStart = rzp[3]
  u = @SVector [rzp[1], rzp[2]]
  phiSpan = (phiStart,phiEnd)
  params = InterpolationParameters(itp, zeros(T,3), 2π/itp.nfp)
  prob = ODEProblem(fieldDerivPhi, u, phiSpan, params)
  if poincare
    N = abs(phiEnd - phiStart)/(2π/itp.nfp)
    saveat = [i * 2*π/itp.nfp + phiStart for i in 1:N]
  else
    saveat = []
  end
  sol = abs(phiStep) > zero(T) ? solve(prob, Tsit5(), dtmax = phiStep, saveat = saveat) :
                                 solve(prob, Tsit5(), saveat = saveat)

end

function followFieldS(itp::BFieldInterpolator{T},
                      rzp::Array{Float64},
                      sEnd::Float64;
                      sStep::Float64=zero(T),
                      ) where{T}
  u = @SVector [rzp[1], rzp[2], rzp[3]]
  sSpan = (0, sEnd)
  params = InterpolationParameters(itp, zeros(T,3), 2π/itp.nfp)
  prob = ODEProblem(fieldDerivS, u, sSpan, params)
  sol = abs(sStep) > zero(T) ? solve(prob, Tsit5(), dtmax = sStep, saveat = sStep) :
                               solve(prob, Tsit5(), saveat = sStep)
end

function fieldDerivPhi!(du::Vector{Float64}, u::Vector{Float64},
                        itp::BFieldInterpolator, ϕ::Float64)
  r = u[1]
  z = u[2]
  br, bz, bp = itp(r, z, ϕ)
  du[1] = r * br/bp
  du[2] = r * bz/bp
end

function fieldDerivPhi(u::AbstractVector,
                       itp::BFieldInterpolator{T},
                       ϕ::Float64;
                      ) where {T}
  ϕ = mod(ϕ, 2π/itp.nfp)
  br, bz, bp = itp(u[1], u[2], ϕ)
  dr = u[1] * br/bp
  dz = u[1] * bz/bp
  SVector{2,T}(dr, dz)
end


function fieldDerivPhi(u::AbstractVector{T},
                       p::InterpolationParameters{T},
                       ϕ::T;
                      ) where {T}
  ϕ = mod(ϕ, p.ϕ_max)
  map!(i->getfield(p.itp, i)(u[1], u[2], ϕ), p.values, 1:3)
  bp = u[1] / p.values[3]
  SVector{2,T}(bp * p.values[1], bp * p.values[2])
end


#integration with respect to arclength
function fieldDerivS(u::AbstractVector,
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
  prob = iszero(ϕ_step) ? ODEProblem(fieldDerivPhi, u0, ϕ_span, params, saveat = ϕ_saveat) :
                          ODEProblem(fieldDerivPhi, u0, ϕ_span, params, saveat = ϕ_saveat, dtmax = ϕ_step)

  function prob_func(prob, i, repeat)
    @info "Remaking at iteration $i with initial condition $(initial_conditions[i][1:2])"
    remake(prob, u0 = SVector{2,T}(initial_conditions[i][1:2]))
  end

  poincare_prob = EnsembleProblem(prob, prob_func = prob_func)
  traj = solve(poincare_prob, Tsit5(), EnsembleThreads(), trajectories = length(initial_conditions))
end
