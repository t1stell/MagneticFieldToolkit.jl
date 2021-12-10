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
    N = abs(phiEnd - phiStart)/(2*π/itp.nfp)
    saveat = [i * 2*π/itp.nfp + phiStart for i in 1:N]
  else
    saveat = []
  end
  sol = abs(phiStep) > zero(T) ? solve(prob, Tsit5(), dtmax = phiStep, saveat = saveat) :
                                 solve(prob, Tsit5(), saveat = saveat)

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

function poincare(itp::BFieldInterpolator,
                  bfield::BField,
                  phi_end::T;
                  phi_step::T=zero(T),
                 ) where {T}
  phi_span = (zero(T), phi_end)
  
end
