


function followField(itp::BFieldInterpolator,
                     rzp::Array{Float64},
                     phiEnd::Float64,
                     phiStep::Float64;
                     poincare::Bool=false
                    )
  phiStart = rzp[3]
  u = @SVector [rzp[1], rzp[2]]
  phiSpan = (0.0,phiEnd)
  prob = ODEProblem(fieldDerivPhi, u, phiSpan, itp)
  sol = solve(prob, Tsit5(), dtmax=phiStep, saveat = poincare ? 2π : [])

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
                       itp::BFieldInterpolator,
                       ϕ::Float64;
                      )
  br, bz, bp = itp(u[1], u[2], ϕ)
  dr = u[1] * br/bp
  dz = u[1] * bz/bp
  @SVector [dr, dz]
end
