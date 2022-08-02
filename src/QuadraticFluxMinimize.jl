"""
  function to calculate the toroidal flux = A ⋅ t
  = dR/dθ * Ar + dZ/dθ * Az
"""
function toroidal_flux(surf::Union{FourierSurface, VmecSurface},
                  field::MagneticField;
                  resolution = 100)
  θ = range(0, 2π, resolution)
  fcs = [FourierCoordinates(1.0, θi, 0.0) for θi in θ]
  R = [inverseTransform(x, surf.rmn) for x in fcs]
  Z = [inverseTransform(x, surf.zmn) for x in fcs]

  BandA = [field(R[i],0.0,Z[i], A=true) for i in 1:length(R)]
  Ar = [BandA[i][2][1] for i in 1:length(BandA)]
  Az = [BandA[i][2][3] for i in 1:length(BandA)]

  dRdθ = [inverseTransform(x, surf.rmn; deriv=:dθ) for x in fcs]
  dZdθ = [inverseTransform(x, surf.zmn; deriv=:dθ) for x in fcs]
  emag = sqrt.(dRdθ.^2 .+ dZdθ.^2)
  tR = dRdθ./emag
  tZ = dZdθ./emag

  #need to integrate here
  dsdθ = [sqrt((R[i+1] - R[i])^2 + (Z[i+1] - Z[i])^2) for i in 1:length(R)-1]
  append!(dsdθ, dsdθ[1])
  return sum(dsdθ.*(Ar.*tR + Az.*tZ))
  


end
