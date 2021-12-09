


function followField(itp::BFieldInterpolator, rzp::Array{Float64}, 
                     phiEnd::Float64, phiStep::Float64)
  phiStart = rzp[3]
  u = [rzp[1], rzp[2]]
  phiRange = phiStart:phiStep:phiEnd
  rvals = Array{Float64}(undef, length(phiRange)+1)
  zvals = Array{Float64}(undef, length(phiRange)+1)
  pvals = Array{Float64}(undef, length(phiRange)+1)
  rvals[1] = rzp[1]
  zvals[1] = rzp[2]
  pvals[1] = rzp[3]
  count = 1
  println("len: ",length(rvals))
  for ϕ in phiRange
    count += 1
    phiSpan = (ϕ, ϕ+phiStep)
    println("at phi: ",ϕ," count ",count, " phiSpan ",phiSpan)
    prob = ODEProblem(fieldDerivPhi!,u,phiSpan,itp)
    sol = solve(prob,Tsit5())

    rvals[count] = sol.u[end][1]
    zvals[count] = sol.u[end][2]
    pvals[count] = sol.t[end]
  end

  return rvals, zvals, pvals



end


function fieldDerivPhi!(du::Array{Float64}, u::Array{Float64}, 
                        itp::BFieldInterpolator, ϕ::Float64)
  r = u[1]
  z = u[2]
  br, bz, bp = itp(r, z, ϕ)
  du[1] = r * br/bp
  du[2] = r * bz/bp
end
