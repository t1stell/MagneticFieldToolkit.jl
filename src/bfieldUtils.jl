"""
  compute_magnetic_potential(cset, xyz; rtol=rtol, currents=currents)
  compute_magnetic_potential(cset, cc; rtol=rtol, currents=currents)
  compute_magnetic_potential(coil, xyz; rtol=rtol, current=current)
  compute_magnetic_potential(coil, cc; rtol=rtol, current=current)

Function to compute a magnetic potential at a point from either an individual coil or a set of coils

# Arguments
 - `cset::Coilset`: A set of coils, can be obtained from using read_vmec_coils on an ASCII coils file
 - `coil::CoilFilamet`: A singular filament (no standalone constructor yet)
 - `xyz::SVector(3)`: A StaticVector representing the point to calculate the potential at in Cartesian coordinates
 - `cc::Cylindrical`: A StaticVector representing the point to calculate the potential at in Cylindrical coordinates

# Optional Arguments
 - `rtol::Float64`: The required tolerance, default is 1.0E-10
 - `currents::Vector{Float64}` or `current::Float64`: for a coil set, this is an array of currents, one for each family.  For a individual coil, it's the current in that coil.

# Outputs
 - `A::SVector(3)`: The vector potential field obtained from integrating μ_0/4π * I * ∫ds⃗/ξ where I is the current, and ξ is the distance between the point on the coil and the target point.

"""
function compute_magnetic_potential(cset::CoilSet{T}, xyz::SVector;
                            rtol=1.0E-10, currents=nothing, exact=false) where {T}
  A = [0.0, 0.0, 0.0]
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  for (idx, family) in enumerate(cset.family)
    if currents == nothing
      current = nothing
    else
      current = currents[idx]
    end
    for coil in family.coil
      At = compute_magnetic_potential(coil, xyz, rtol=rtol, current=current, exact=exact)
      A = A .+ At
    end
  end
  return SVector(A)
end

function compute_magnetic_potential(cset::CoilSet{T}, cc::Cylindrical;
                            rtol = 1.0E-10, currents=nothing,
                            exact=false) where {T}
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  xyz = CartesianFromCylindrical()(cc)
  (Ax, Ay, Az) =  compute_magnetic_potential(cset, xyz, rtol=rtol, currents=currents, exact=exact)
  #convert to Cylindrical
  Ar = Ax * cos(cc.θ) + Ay * sin(cc.θ)
  Aθ = -Ax * sin(cc.θ) + Ay * cos(cc.θ)
  return SVector(Ar, Aθ, Az)
end

function compute_magnetic_potential(coil::CoilFilament{T}, xyz::SVector;
                            rtol = 1.0E-10, current=nothing,
                            exact=false) where{T}
  if exact
    return compute_magnetic_potential_exact(coil, xyz, rtol=rtol, current=current)
  end
  xch = coil.xnodes_half
  ych = coil.ynodes_half
  zch = coil.znodes_half
  dxc = coil.dxnodes
  dyc = coil.dynodes
  dzc = coil.dznodes
  seg_res = length(xch)
  ξ_x = [xyz[1] - xch[i] for i in 1:seg_res]
  ξ_y = [xyz[2] - ych[i] for i in 1:seg_res]
  ξ_z = [xyz[3] - zch[i] for i in 1:seg_res]
  ξ_mag = sqrt.(ξ_x.^2 .+ ξ_y.^2 .+ ξ_z.^2)
  A = sum([ SVector(dxc[i], dyc[i], dzc[i]) / ξ_mag[i] for i in 1:seg_res])
  if current == nothing
    current = coil.current
  end
  return A * current * μ0over4π 
  
end



function compute_magnetic_potential(coil::CoilFilament{T}, cc::Cylindrical;
                            rtol = 1.0E-10, current=nothing, exact=false) where {T}
  xyz = CartesianFromCylindrical()(cc)
  (Ax, Ay, Az) =  compute_magnetic_potential(coil, xyz, rtol=rtol, current=current, exact=exact)
  #convert to Cylindrical
  Ar = Ax * cos(cc.θ) + Ay * sin(cc.θ)
  Aθ = -Ax * sin(cc.θ) + Ay * cos(cc.θ)
  return SVector(Ar, Aθ, Az)
end

function compute_magnetic_potential_exact(coil::CoilFilament{T}, xyz::SVector;
                                          rtol = 1.0E-10, current=nothing) where{T}
  xc = coil.x
  yc = coil.y
  zc = coil.z
  function get_A(t::Float64)
    ξ = sqrt((xyz[1] - xc(t))^2 + (xyz[2] -yc(t))^2 + (xyz[3] - zc(t))^2)
    return SVector(coil.dxdt(t)/(coil.ds_mag(t) * ξ),
                   coil.dydt(t)/(coil.ds_mag(t) * ξ),
                   coil.dzdt(t)/(coil.ds_mag(t) * ξ))
  end
  if current == nothing
    current = coil.current
  end
  A = hquadrature(get_A, 0, 2π, rtol=1.0E-10)[1]
  return SVector(A .* current * μ0over4π)
end

"""
  compute_magnetic_field(cset, xyz; rtol=rtol, currents=currents)
  compute_magnetic_field(cset, cc; rtol=rtol, currents=currents)
  compute_magnetic_field(coil, xyz; rtol=rtol, current=current)
  compute_magnetic_field(coil, cc; rtol=rtol, current=current)

Function to compute a magnetic field at a point from either an individual coil or a set of coils

# Arguments
 - `cset::Coilset`: A set of coils, can be obtained from using read_vmec_coils on an ASCII coils file
 - `coil::CoilFilamet`: A singular filament (no standalone constructor yet)
 - `xyz::SVector(3)`: A StaticVector representing the point to calculate the field at in Cartesian coordinates
 - `cc::Cylindrical`: A StaticVector representing the point to calculate the field at in Cylindrical coordinates

# Optional Arguments
 - `rtol::Float64`: The required tolerance, default is 1.0E-10
 - `currents::Vector{Float64}` or `current::Float64`: for a coil set, this is an array of currents, one for each family.  For a individual coil, it's the current in that coil.

# Outputs
 - `B::SVector(3)`: The vector magnetic field obtained from integrating μ_0/4π * I * ∫ds⃗ x ξ⃗/ξ^2 where I is the current, and ξ is the vector between the point on the coil and the target point.

"""
function compute_magnetic_field(cset::CoilSet{T}, xyz::SVector;
                            rtol=1.0E-10, currents=nothing, exact=false) where {T}
  B = [0.0, 0.0, 0.0]
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  for (idx, family) in enumerate(cset.family)
    if currents == nothing
      current = nothing
    else
      current = currents[idx]
    end
    for coil in family.coil
      Bt = compute_magnetic_field(coil, xyz, rtol=rtol, current=current, exact=exact)
      B = B .+ Bt
    end
  end
  return SVector(B)
end

function compute_magnetic_field(cset::CoilSet{T}, cc::Cylindrical;
                         rtol=1.0E-10, currents=nothing, exact=false) where{T}
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  xyz = CartesianFromCylindrical()(cc)
  (Bx, By, Bz) =  compute_magnetic_field(cset, xyz, rtol=rtol, currents=currents, exact=exact)
  #convert to Cylindrical
  Br = Bx * cos(cc.θ) + By * sin(cc.θ)
  Bθ = -Bx * sin(cc.θ) + By * cos(cc.θ)
  return SVector(Br, Bθ, Bz)
end

function compute_magnetic_field(cset::CoilSet{T}, coords::C;
                                rtol = 1.0E-10, currents=nothing, exact=false, xyz=false
                               ) where {T, C<:AbstractArray}
    if xyz
        xyzc = SVector{3}(coords)
        return compute_magnetic_field(cset, xyzc, rtol=rtol, currents=currents, exact=exact)
    else #in rθz
        rθz = Cylindrical(coords[1], coords[2], coords[3])
        return compute_magnetic_field(cset, rθz, rtol=rtol, currents=currents, exact=exact)
    end
end


function compute_magnetic_field(coil::CoilFilament{T}, xyz::SVector;
                         rtol=1.0E-10, current=nothing, 
                         exact=false) where{T}
  if exact
    return compute_magnetic_field_exact(coil, xyz, rtol=rtol, current=current)
  end
  xch = coil.xnodes_half
  ych = coil.ynodes_half
  zch = coil.znodes_half
  dxc = coil.dxnodes
  dyc = coil.dynodes
  dzc = coil.dznodes
  seg_res = length(xch)
  ξ_x = [xyz[1] - xch[i] for i in 1:seg_res]
  ξ_y = [xyz[2] - ych[i] for i in 1:seg_res]
  ξ_z = [xyz[3] - zch[i] for i in 1:seg_res]
  ξ_cubed = (ξ_x.^2 .+ ξ_y.^2 .+ ξ_z.^2).^(1.5)
  dlcrossr = [cross(SVector(dxc[i], dyc[i], dzc[i]), SVector(ξ_x[i], ξ_y[i], ξ_z[i])) 
             for i in 1:seg_res]
  if current == nothing
    current = coil.current
  end
  return (current * μ0over4π)*sum(dlcrossr ./ ξ_cubed)

end

function compute_magnetic_field(coil::CoilFilament{T}, cc::Cylindrical;
                         rtol=1.0E-10, current=nothing, exact=false) where{T}
  xyz = CartesianFromCylindrical()(cc)
  (Bx, By, Bz) =  compute_magnetic_field(coil, xyz, rtol=rtol, current=current, exact=exact)
  #convert to Cylindrical
  Br = Bx * cos(cc.θ) + By * sin(cc.θ)
  Bθ = -Bx * sin(cc.θ) + By * cos(cc.θ)
  return SVector(Br, Bθ, Bz)
end

function compute_magnetic_field(coil::CoilFilament{T}, coords::C;
                         rtol = 1.0E-10, current=nothing, exact=false, xyz=false
                               ) where {T, C <: AbstractArray}
    if xyz
        xyzc = SVector{3}(coords)
        return compute_magnetic_field(coil, xyzc, rtol=rtol, current=current, exact=exact)
    else #in rθz
        rθz = Cylindrical(coords[1], coords[2], coords[3])
        return compute_magnetic_field(coil, rθz, rtol=rtol, current=current, exact=exact)
    end
end

#More accurate integration using quadrature.  Slow
function compute_magnetic_field_exact(coil::CoilFilament{T}, xyz::SVector;
                                      rtol=1.0E-10, current=nothing) where {T}
  xc = coil.x
  yc = coil.y
  zc = coil.z
  function get_B(t::Float64)
    
    ξ = SVector((xyz[1] - xc(t)),(xyz[2] -yc(t)),(xyz[3] - zc(t)))
    ξ_cubed = (ξ[1]^2 + ξ[2]^2 + ξ[3]^2).^(1.5)
    dlcrossr = cross(SVector(coil.dxdt(t), coil.dydt(t), coil.dzdt(t)), ξ)
    return SVector(dlcrossr ./ ξ_cubed)
  end
  if current == nothing
    current = coil.current
  end
  B = hquadrature(get_B, 0, 2π, rtol=rtol)[1]
  return SVector(B .* (current * μ0over4π))

end
