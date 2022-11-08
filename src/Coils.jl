const μ0over4π = 1.25663706E-6/4/π

struct CoilFilament{FT}
  x::Interpolations.Extrapolation
  y::Interpolations.Extrapolation
  z::Interpolations.Extrapolation
  r::Interpolations.Extrapolation
  dxdt::Interpolations.Extrapolation
  dydt::Interpolations.Extrapolation
  dzdt::Interpolations.Extrapolation
  drdt::Interpolations.Extrapolation
  ds_mag::Interpolations.Extrapolation
  current::FT
end

struct CoilFamily{FT}
  coil::Array{CoilFilament{FT}, 1}
  name::String
end
"""
  The coil hierarchy uses the convention common in various coil formats, like
  the ones used by vmec and mgrid. 
  Each coil set is broken into a family of coil filaments that all have
  the same currents, and each family is broken up into individual filaments
  Accessing a specific coil is done by:
  coil_set.family[i].coil[j]. etc
"""
struct CoilSet{FT}
  family::Array{CoilFamily{FT},1}
end

"""
  read_vmec_coils(filename, use_current=false)

Reads a coil file in the vmec format and outputs a coil_set 
Right now the code assumes all coils exist and no explicit rotations
are calculated. Check coil file to verify.

Coil files tend to have assigned currents. Typically we want to set these
all to zero, so we can specify our own currents. If you want to use the
actual currents in the file (equivalent to the "raw current") option in
mgrid, you can set the "use_current" flag to true.
"""
function read_vmec_coils(filename::String; use_current = false)
  #we read through the file twice, first to determine how many coil families
  #there are and how many coils in each family, then again
  #to save the actual coil data
  family_names = Vector{String}(undef, 0)
  coils_per_family = Vector{Int64}(undef, 0)
  lines = readlines(filename)
  coilends = Vector{Int64}(undef, 0)
  family_labels = Vector{Int64}(undef, 0)
  for (i, ln) in enumerate(lines)
    dum = split(ln)
    if length(dum) < 5
      continue
    end
    family_name = dum[6]
    family_label = parse(Int, dum[5])
    while family_label > length(coils_per_family)
      push!(coils_per_family, 0)
      push!(family_names, "")
    end
    coils_per_family[family_label] += 1
    if family_names[family_label] == ""
      family_names[family_label] = family_name
    end
    push!(coilends, i)
    push!(family_labels, family_label) 
  end
  nfamilies = length(family_names)
  nfilaments = length(coilends)
  
  #Now read back through the file and assign stuff
  coil_set = Vector{CoilFamily{Float64}}(undef, nfamilies)
  for fam_index in 1:nfamilies
    coil_family = Vector{CoilFilament{Float64}}(undef, coils_per_family[fam_index])
    #get all the indices
    family_indices = findall(x->x==fam_index, family_labels)
    for (j, coil_index) in enumerate(family_indices)
      #guess at start index
      coil_index == 1 ? start_index = 1 : start_index = coilends[coil_index-1]+1
      end_index = coilends[coil_index]
      #we can't preassign the vectors bc there might be junk values
      xc = Vector{Float64}(undef, 0)
      yc = Vector{Float64}(undef, 0)
      zc = Vector{Float64}(undef, 0)
      current = nothing
      for j in start_index:end_index
        dum = split(lines[j])
        if length(dum) < 4
          continue
        end
        push!(xc, parse(Float64, dum[1]))
        push!(yc, parse(Float64, dum[2]))
        push!(zc, parse(Float64, dum[3]))
        if current == nothing && use_current
          current = parse(Float64, dum[4])
        else
          current = 1.0 #scaled current
        end
      end
      #make the splines
      ts = range(0, 2π, length(xc))
      rc = sqrt.(xc.^2 + yc.^2)
      xs = cubic_spline_interpolation(ts, xc)
      ys = cubic_spline_interpolation(ts, yc)
      zs = cubic_spline_interpolation(ts, zc)
      rs = cubic_spline_interpolation(ts, rc)

      #calculate the derivative splines
      dx = [Interpolations.gradient(xs, t)[1] for t in ts] 
      dy = [Interpolations.gradient(ys, t)[1] for t in ts]
      dz = [Interpolations.gradient(zs, t)[1] for t in ts]
      dr = [Interpolations.gradient(rs, t)[1] for t in ts]
      ds_mag = sqrt.(dx.^2 .+ dy.^2 .+ dz.^2)
      dxdt = cubic_spline_interpolation(ts, dx)
      dydt = cubic_spline_interpolation(ts, dy)
      dzdt = cubic_spline_interpolation(ts, dz)
      drdt = cubic_spline_interpolation(ts, dr)
      ds = cubic_spline_interpolation(ts, ds_mag)
      coil_family[j] = CoilFilament(xs, ys, zs, rs, dxdt, dydt, 
                                    dzdt, drdt, ds, current)
    end
    coil_set[fam_index] = CoilFamily(coil_family, family_names[fam_index])
  end
  return CoilSet(coil_set)
  
end 

#function to get approximate extremes of coils to bound box fields
function extreme_coils(cset::CoilSet{T}, coord::Symbol; vmax=true) where T
  if vmax
    vextreme = -1.0E30
  else
    vextreme = 1.0E30 #just pick an impossibly large number
  end
  for family in cset.family
    for coil in family.coil
      vs = getfield(coil, coord)
      v = vs.itp.itp.coefs[1:end-1]
      if vmax
        vextreme_temp = maximum(v)
        if vextreme_temp > vextreme
          vextreme = vextreme_temp
        end
      else
        vextreme_temp = minimum(v)
        if vextreme_temp < vextreme
          vextreme = vextreme_temp
        end
      end
    end
  end
  return vextreme
end

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
                            rtol=1.0E-10, currents=nothing) where {T}
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
      At = compute_magnetic_potential(coil, xyz, rtol=rtol, current=current)
      A = A .+ At
    end
  end
  return SVector(A)
end

function compute_magnetic_potential(cset::CoilSet{T}, cc::Cylindrical;
                            rtol = 1.0E-10, currents=nothing) where {T}
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  xyz = CartesianFromCylindrical()(cc)
  (Ax, Ay, Az) =  compute_magnetic_potential(cset, xyz, rtol=rtol, currents=currents)
  #convert to Cylindrical
  Ar = Ax * cos(cc.θ) + Ay * sin(cc.θ)
  Aθ = -Ax * sin(cc.θ) + Ay * cos(cc.θ)
  return SVector(Ar, Aθ, Az)
end

function compute_magnetic_potential(coil::CoilFilament{T}, xyz::SVector;
                            rtol = 1.0E-10, current=current) where{T}
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

function compute_magnetic_potential(coil::CoilFilament{T}, cc::Cylindrical;
                            rtol = 1.0E-10, current = current) where {T}
  xyz = CartesianFromCylindrical()(cc)
  (Ax, Ay, Az) =  compute_magnetic_potential(coil, xyz, rtol=rtol, current=current)
  #convert to Cylindrical
  Ar = Ax * cos(cc.θ) + Ay * sin(cc.θ)
  Aθ = -Ax * sin(cc.θ) + Ay * cos(cc.θ)
  return SVector(Ar, Aθ, Az)
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
                            rtol=1.0E-10, currents=nothing) where {T}
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
      Bt = compute_magnetic_field(coil, xyz, rtol=rtol, current=current)
      B = B .+ Bt
    end
  end
  return SVector(B)
end

function compute_magnetic_field(cset::CoilSet{T}, cc::Cylindrical;
                         rtol=1.0E-10, currents=nothing) where{T}
  if currents != nothing
    if length(currents) < length(cset.family)
      print("Currents vector must have at least as many entries as coil families")
      return 0.0
    end
  end
  xyz = CartesianFromCylindrical()(cc)
  (Bx, By, Bz) =  compute_magnetic_field(cset, xyz, rtol=rtol, currents=currents)
  #convert to Cylindrical
  Br = Bx * cos(cc.θ) + By * sin(cc.θ)
  Bθ = -Bx * sin(cc.θ) + By * cos(cc.θ)
  return SVector(Br, Bθ, Bz)
end

function compute_magnetic_field(coil::CoilFilament{T}, xyz::SVector;
                         rtol=1.0E-10, current=nothing) where{T}
  xc = coil.x
  yc = coil.y
  zc = coil.z
  function get_B(t::Float64)
    
    ξ = SVector((xyz[1] - xc(t)),(xyz[2] -yc(t)),(xyz[3] - zc(t)))
    ξ_squared = ξ[1]^2 + ξ[2]^2 + ξ[3]^2
    dlcrossr = cross(SVector(coil.dxdt(t), coil.dydt(t), coil.dzdt(t)), ξ)
    return SVector(dlcrossr ./ ξ_squared)
  end
  if current == nothing
    current = coil.current
  end
  B = hquadrature(get_B, 0, 2π, rtol=rtol)[1]
  return SVector(B .* (current * μ0over4π))

end

function compute_magnetic_field(coil::CoilFilament{T}, cc::Cylindrical;
                         rtol=1.0E-10, current=nothing) where{T}
  xyz = CartesianFromCylindrical()(cc)
  (Bx, By, Bz) =  compute_magnetic_field(coil, xyz, rtol=rtol, current=current)
  #convert to Cylindrical
  Br = Bx * cos(cc.θ) + By * sin(cc.θ)
  Bθ = -Bx * sin(cc.θ) + By * cos(cc.θ)
  return SVector(Br, Bθ, Bz)
end


"""
  potential_from_coils(cset::Coilset, r_res::Int, z_res::Int)

optional arguments, rmin, rmax, zmin, zmax.  If not set, the values used are
the maximum r and z from the coils

"""
function potential_from_coils(cset::CoilSet{T}, r_res::Int64, z_res::Int64,
                               nfp::Int64;
                  rmin=nothing, rmax=nothing, zmin=nothing, zmax=nothing
                  ) where {T}
  if rmax == nothing                
    rmax = extreme_coils(cset, :r, vmax=true)
  end
  if rmin == nothing
    rmin = extreme_coils(cset, :r, vmax=false)
  end
  if zmax == nothing
    zmax = extreme_coils(cset, :z, vmax = true)
  end
  if zmin == nothing
    zmin = extreme_coils(cset, :z, vmin = true)
  end

  #calculate the potential here
end  



