using HDF5
const μ0over4π = 1.25663706E-6/4/π


#this is an mgrid that allows for multiple values of the field
#so currents can be adjusted independently
struct MultipleFieldGrid <: AbstractMagneticField
  magnetic_field::Vector{MagneticField}
end

#constructor to generate a field grid
function MultipleFieldGrid(r::StepRangeLen, θ::StepRangeLen, 
                           z::StepRangeLen, nfp::Int, nfamilies::Int) 
  dummy = zeros(length(r), length(θ), length(z))

  #construct the coords
  fullSize = (length(r), length(θ), length(z))
  r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
  θ_grid = reshape(repeat(θ,inner=length(r),outer=length(z)),fullSize)
  z_grid = reshape(repeat(z,inner=length(r)*length(θ)),fullSize)
  field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))

  
  gridarray = fill(MagneticField(field_coords, copy(dummy), copy(dummy), copy(dummy),
                            copy(dummy), copy(dummy), copy(dummy), nfp=nfp), nfamilies)

  return MultipleFieldGrid(gridarray)
end

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
  #maybe change these to SVectors?
  xnodes::Array{FT, 1}
  ynodes::Array{FT, 1}
  znodes::Array{FT, 1}
  rnodes::Array{FT, 1}
  #compute on half nodes for quick biot-savart calculation
  xnodes_half::Array{FT, 1}
  ynodes_half::Array{FT, 1}
  znodes_half::Array{FT, 1}
  rnodes_half::Array{FT, 1}
  dxnodes::Array{FT, 1}
  dynodes::Array{FT, 1}
  dznodes::Array{FT, 1}
  drnodes::Array{FT, 1}
  dsnodes::Array{FT, 1}
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

function generate_coil_filament(xc::T, yc::T, zc::T, current::F;
                                node_res=1025) where {T, F}
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
  
  #calculate coils on the node resolution
  ts = range(0, 2π, node_res)
  xnodes = [xs(t) for t in ts]
  ynodes = [ys(t) for t in ts]
  znodes = [zs(t) for t in ts]
  rnodes = [rs(t) for t in ts]
  #compute the half nodes
  tstep = ts[2] - ts[1]
  tshalf = range(tstep/2, 2π - tstep/2, node_res-1)
  xnodes_half = [xs(t) for t in tshalf]
  ynodes_half = [ys(t) for t in tshalf]
  znodes_half = [zs(t) for t in tshalf]
  rnodes_half = [rs(t) for t in tshalf]
  # compute the lengths along the coils
  dsnodes = [ sqrt( (xnodes[i+1] - xnodes[i])^2 + 
                    (ynodes[i+1] - ynodes[i])^2 + 
                    (znodes[i+1] - znodes[i])^2) for i in 1:node_res-1]
  dxnodes = [(xnodes[i+1] - xnodes[i]) for i in 1:node_res-1]
  dynodes = [(ynodes[i+1] - ynodes[i]) for i in 1:node_res-1]
  dznodes = [(znodes[i+1] - znodes[i]) for i in 1:node_res-1]
  drnodes = [(rnodes[i+1] - rnodes[i]) for i in 1:node_res-1]

  return CoilFilament(xs, ys, zs, rs, dxdt, dydt, 
                      dzdt, drdt, ds, 
                      xnodes, ynodes, znodes, rnodes,
                      xnodes_half, ynodes_half, znodes_half, rnodes_half,
                      dxnodes, dynodes, dznodes, drnodes,
                      dsnodes, current)
end

function read_coil_files(filename_list::Vector{String}, nfp::Integer;
                         grouping::Vector{Int64}=[1], currents::Vector{Float64}=[1.0])
  nfamilies = maximum(grouping)
  nunique = length(filename_list)
  nfilaments = nunique*nfp
  coil_set = Vector{CoilFamily{Float64}}(undef, nfamilies)
  #for now assume just one family
  coil_family = Vector{CoilFilament{Float64}}(undef, nfilaments)
  for (file_ix, filename) in enumerate(filename_list)
    lines = readlines(filename)
    npoints = parse(Int64, lines[1])
    xc = Vector{Float64}(undef, npoints)
    yc = Vector{Float64}(undef, npoints)
    zc = Vector{Float64}(undef, npoints)
    for (i, line) in enumerate(lines[2:npoints+1])
      dum = split(line)
      xc[i] = parse(Float64, dum[1])
      yc[i] = parse(Float64, dum[2])
      zc[i] = parse(Float64, dum[3])
    end
    coil_family[file_ix] = generate_coil_filament(xc, yc, zc, currents[1])
    #now do the rotations

    rc = sqrt.(xc.^2 + yc.^2)
    θ = atan.(yc, xc)
    for n in 2:nfp
      θshift = (2π / nfp) * (n-1)
      θn = θ .+ θshift
      xn = rc .* cos.(θn)
      yn = rc .* sin.(θn)
      coil_family[file_ix + (nunique * (n-1))] = generate_coil_filament(xn, yn, zc, currents[1])
    end
  end
  return CoilSet([CoilFamily(coil_family, "coil")])
end

"""
  read_vmec_coils(filename, use_current=false)

Reads a coil file in the vmec format and outputs a coil_set 
Right now the code assumes all coils exist and no explicit rotations
are calculated. Check coil file to verify.

Coil files tend to have assigned currents. Typically we want to set these
all to one, so we can specify our own currents. If you want to use the
actual currents in the file (equivalent to the "raw current") option in
mgrid, you can set the "use_current" flag to true.
"""
function read_vmec_coils(filename::String; use_current = false, node_res = 1025)
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
    for (k, coil_index) in enumerate(family_indices)
      #guess at start index
      coil_index == 1 ? start_index = 1 : start_index = coilends[coil_index-1]+1
      end_index = coilends[coil_index]
      #we can't preassign the vectors bc there might be junk values
      xc = Vector{Float64}(undef, 0)
      yc = Vector{Float64}(undef, 0)
      zc = Vector{Float64}(undef, 0)
      
      current = nothing
      shouldflip = false
      for j in start_index:end_index
        dum = split(lines[j])
        if length(dum) < 4
          continue
        end
        push!(xc, parse(Float64, dum[1]))
        push!(yc, parse(Float64, dum[2]))
        push!(zc, parse(Float64, dum[3]))
        tempcurrent = parse(Float64, dum[4])
        if j == end_index-1 && tempcurrent < 0 
          shouldflip = true
        end
        if j == end_index-1 && use_current
          current = tempcurrent
        else
          current = 1.0 #scaled current
        end
      end
      #flip the reversed ones (this is a really dumb convention)
      if shouldflip
        xc = reverse(xc)
        yc = reverse(yc)
        zc = reverse(zc)
      end

      coil_family[k] = generate_coil_filament(xc, yc, zc, current, node_res=node_res)
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

#helper function for generate mgrid
function get_bfield_grid!(coil::CoilFilament{T},
                         Br::Array{T}, Bz::Array{T}, 
                         Bϕ::Array{T}, Ar::Array{T},
                         Az::Array{T}, Aϕ::Array{T},
                         rrange::StepRangeLen, zrange::StepRangeLen, 
                         ϕrange::StepRangeLen) where {T}
  z_res = length(zrange)
  ϕ_res = length(ϕrange)

  #handle difference in stell symmetry between odd/even phi res
  if ϕ_res % 2 == 0
    ϕ_half = div(ϕ_res, 2)
    ϕodd = false
  else
    ϕ_half = div(ϕ_res,2) + 1
    ϕodd = true
  end



  for (ϕ_i) in 1:ϕ_half #only need to do the half period
    ϕ = ϕrange[ϕ_i]
    ϕ_is = ϕ_res - ϕ_i + 1 # stell symmetric ϕ 
    for (z_i,z) in enumerate(zrange)
      z_is = z_res - z_i + 1 # stell symmetric z
      for (r_i, r) in enumerate(rrange)
        cc = Cylindrical(r, ϕ, z)
        Bpoint = compute_magnetic_field(coil, cc, current=1.0)
        Apoint = compute_magnetic_potential(coil, cc, current=1.0)
        #assign the grid points
        #use the mgrid convention for ordering
        Br[r_i, ϕ_i, z_i] += Bpoint[1]
        Bz[r_i, ϕ_i, z_i] += Bpoint[3]
        Bϕ[r_i, ϕ_i, z_i] += Bpoint[2]
        Ar[r_i, ϕ_i, z_i] += Apoint[1]
        Aϕ[r_i, ϕ_i, z_i] += Apoint[3]
        Az[r_i, ϕ_i, z_i] += Apoint[2]
        #don't do stell sym for this surface
        if ϕodd && ϕ_i == ϕ_half
          continue
        end
        Br[r_i, ϕ_is, z_is] += -Bpoint[1]
        Bz[r_i, ϕ_is, z_is] += Bpoint[3]
        Bϕ[r_i, ϕ_is, z_is] += Bpoint[2]
        Ar[r_i, ϕ_is, z_is] += -Apoint[1]
        Az[r_i, ϕ_is, z_is] += Apoint[3]
        Aϕ[r_i, ϕ_is, z_is] += Apoint[2]
      end
    end
  end
end

"""
  generate_mgrid(cset::Coilset, r_res::Int, z_res::Int, ϕ_res::Int, nfp::Int, savename)

optional arguments, rmin, rmax, zmin, zmax.  If not set, the values used are
the maximum r and z from the coils

The code assumes stellarator symmetry

"""
function generate_mgrid(cset::CoilSet{T}, r_res::Int64, θ_res::Int64, 
          z_res::Int64,  nfp::Int64; 
                  savename = nothing,
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
    zmin = extreme_coils(cset, :z, vmax = false)
  end

  if zmax != -1 * zmin
    #stell symmetry won't work if zmin and zmax aren't the same
    zmin = -1*zmax
  end

  Br_temp = Array{Float64}(undef, r_res, θ_res, z_res)  
  Bz_temp = Array{Float64}(undef, r_res, θ_res, z_res)  
  Bθ_temp = Array{Float64}(undef, r_res, θ_res, z_res)  
  Ar_temp = Array{Float64}(undef, r_res, θ_res, z_res)  
  Az_temp = Array{Float64}(undef, r_res, θ_res, z_res)  
  Aθ_temp = Array{Float64}(undef, r_res, θ_res, z_res)  

  #calculate the ranges
  rrange = range(rmin, rmax, r_res)
  zrange = range(zmin, zmax, z_res)
  #note mgrids do not include the last ϕ value because it should
  #be the same, we could do the same thing if we wanted
  θrange = range(0, 2*π/nfp, θ_res)

  #Initialize a magnetic field grid
  mgrid = MultipleFieldGrid(rrange, θrange, zrange, nfp, length(cset.family))
  #generate coords
  fullSize = (length(rrange), length(θrange), length(zrange))
  r_grid = reshape(repeat(rrange,outer=length(zrange)*length(θrange)),fullSize)
  θ_grid = reshape(repeat(θrange,inner=length(rrange),outer=length(zrange)),fullSize)
  z_grid = reshape(repeat(zrange,inner=length(rrange)*length(θrange)),fullSize)
  field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))


  
  # we need to compute for each family
  for (family_idx, family) in enumerate(cset.family)
    #set everything to 0
    Br_temp .= 0.0
    Ar_temp .= 0.0
    Bz_temp .= 0.0
    Az_temp .= 0.0
    Bθ_temp .= 0.0
    Aθ_temp .= 0.0
    count = 0
    for coil in family.coil
      count += 1
      get_bfield_grid!(coil, Br_temp, Bz_temp, Bθ_temp,
                       Ar_temp, Az_temp, Aθ_temp,
                       rrange, zrange, θrange)
      #generate a magnetic field object
      mf = MagneticField(field_coords, Br_temp, Bθ_temp, Bz_temp, Ar_temp, Aθ_temp, Az_temp, 
                         nfp = nfp)
      #load it into the array
      mgrid.magnetic_field[family_idx] = mf
    end
  end
  if savename != nothing
    println("saving the grid, this may take a few minutes")
    save_mgrid(mgrid, savename)
  end

  return mgrid 
end  

function save_mgrid(mgrid::MultipleFieldGrid, savename::String)
  
  fid = h5open(savename, "w")
  k = knots(mgrid.magnetic_field[1].field_data[1].itp)
  r = k.iterators[1].knots
  θ = k.iterators[2].knots
  z = k.iterators[3].knots
  mgsize = size(k)
  ftemp = Array{Float64}(undef, mgsize) #temporary data to store values
  gname = "axes"
  create_group(fid, gname)
  gid = fid[gname]
  write(gid,"r",r)
  write(gid,"theta",θ)
  write(gid,"z",z)
  write(gid,"nfp",mgrid.magnetic_field[1].nfp)
  
  #field arrays
  bnames = ("Br","Btheta","Bz")
  anames = ("Ar","Atheta","Az")

  for ifield in 1:length(mgrid.magnetic_field) #for each field
    gname = "coil"*string(ifield)
    create_group(fid, gname)
    gid = fid[gname]
    #construct the array for magnetic field
    for itype in 1:3
      for ir in 1:mgsize[1]
        for iθ in 1:mgsize[2]
          for iz in 1:mgsize[3]
            ftemp[ir, iθ, iz] = mgrid.magnetic_field[ifield].field_data[itype].itp[ir, iθ, iz]
          end
        end
      end
      write(gid,bnames[itype],ftemp)
    end
    #now do the same for the potential
    for itype in 1:3
      for ir in 1:mgsize[1]
        for iθ in 1:mgsize[2]
          for iz in 1:mgsize[3]
            ftemp[ir, iθ, iz] = mgrid.magnetic_field[ifield].potential_data[itype].itp[ir, iθ, iz]
          end
        end
      end
      write(gid,anames[itype],ftemp)
    end
  end
  close(fid)
end   

   

"""
  functions to call the multiple magnetic field
"""
function (multi_magnetic_field::MultipleFieldGrid)(r::T, θ::T, z::T, currents::Vector{T}) where {T}
  nfields = length(multi_magnetic_field.magnetic_field)
  if nfields != length(currents)
    println("Error: no. of currents, ",length(currents)," does not match no. of fields, ",nfields)
    return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
  end
  Br = 0.0
  Bz = 0.0
  Bθ = 0.0
  Ar = 0.0
  Az = 0.0
  Aθ = 0.0
  for (i, mg) in enumerate(multi_magnetic_field.magnetic_field)
    ((Brt, Bθt, Bzt), (Art, Aθt, Azt)) = mg(r, θ, z, A=true)
    Br += Brt * currents[i]
    Bz += Bzt * currents[i]
    Bθ += Bθt * currents[i]
    Ar += Art * currents[i]
    Az += Azt * currents[i]
    Aθ += Aθt * currents[i]
  end
  return (Br, Bθ, Bz), (Ar, Aθ, Az)
end
