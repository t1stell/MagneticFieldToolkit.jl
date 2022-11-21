
#read an mgrid file and create a relevant struct

#function readMgrid(filename::String)
#  mgrid = NetCDF.open(mgridfile)
#

#read a bmw file and create a relevant field struct
"""
    read_bmw(filename; n_modes=0, padding=5)

Read a NetCDF BMW file and store the result in a `BField` type.

# Arguemnts:
- `filename::AbstractString`   # NetCDF BMW file name

# See also: [`BField`](@ref)
"""
function read_bmw(filename::AbstractString;
                 )
    bmwnetcdf = NetCDF.open(filename)
    bmw_vars = bmwnetcdf.vars
    bmw_dims = bmwnetcdf.dim
    #
    #
    #there's probably a better way to do this, but damn the docs are awful
    nr = 0
    nz = 0
    nθ = 0
    for i in 1:length(bmw_dims.vals)
        #why the fuck there are undefined values here, who knows???
        if isassigned(bmw_dims.vals, i)
            val = bmw_dims.vals[i]
            if val.name == "r"
                nr = Int(val.dimlen)
                @debug "Size of grid in r direction: $nr"
            elseif val.name == "z"
                nz = Int(val.dimlen)
                @debug "Size of grid in z direction: $nz"
            elseif val.name == "phi"
                nθ = Int(val.dimlen)
                @debug "Size of grid in ϕ direction: $nphi"
            end
        end
    end

    rmin = NetCDF.readvar(bmw_vars["rmin"])[1]
    rmax = NetCDF.readvar(bmw_vars["rmax"])[1]
    zmin = NetCDF.readvar(bmw_vars["zmin"])[1]
    zmax = NetCDF.readvar(bmw_vars["zmax"])[1]
    nfp = NetCDF.readvar(bmw_vars["nfp"])[1]
    
    r = range(rmin, rmax, nr)
    z = range(zmin, zmax, nz)
    θ = range(0, 2*π/nfp, nθ+1)

    fullSize = (length(r), length(θ), length(z))
    r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
    θ_grid = reshape(repeat(θ,inner=length(r),outer=length(z)),fullSize)
    z_grid = reshape(repeat(z,inner=length(r)*length(θ)),fullSize)
    field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))

    Br = permutedims(NetCDF.readvar(bmw_vars["br_grid"]),[1,3,2])
    Bz = permutedims(NetCDF.readvar(bmw_vars["bz_grid"]),[1,3,2])
    Bθ = permutedims(NetCDF.readvar(bmw_vars["bp_grid"]),[1,3,2])
    
    Brf = zeros(nr, nθ+1, nz)
    Bzf = zeros(nr, nθ+1, nz)
    Bθf = zeros(nr, nθ+1, nz)

    Brf[:,1:end-1,:] = Br[:,:,:]
    Bzf[:,1:end-1,:] = Bz[:,:,:]
    Bθf[:,1:end-1,:] = Bθ[:,:,:]
    
    Brf[:,end,:] = Br[:,1,:]
    Bzf[:,end,:] = Bz[:,1,:]
    Bθf[:,end,:] = Bθ[:,1,:]


    Ar = permutedims(NetCDF.readvar(bmw_vars["ar_grid"]),[1,3,2])
    Az = permutedims(NetCDF.readvar(bmw_vars["az_grid"]),[1,3,2])
    Aθ = permutedims(NetCDF.readvar(bmw_vars["ap_grid"]),[1,3,2])
    
    Arf = zeros(nr, nθ+1, nz)
    Azf = zeros(nr, nθ+1, nz)
    Aθf = zeros(nr, nθ+1, nz)

    Arf[:,1:end-1,:] = Ar[:,:,:]
    Azf[:,1:end-1,:] = Az[:,:,:]
    Aθf[:,1:end-1,:] = Aθ[:,:,:]
    
    Arf[:,end,:] = Ar[:,1,:]
    Azf[:,end,:] = Az[:,1,:]
    Aθf[:,end,:] = Aθ[:,1,:]

    
    return MagneticField(field_coords, Brf, Bθf, Bzf, Arf, Aθf, Azf, nfp=nfp)

end

# todo: This can be updated to return a multiplegrid object, perhaps with the option
# to condense into one.
function read_mgrid(filename::AbstractString,
                   current::Vector{Float64};
                  )

    mgridnetcdf = NetCDF.open(filename)
    mgrid_vars = mgridnetcdf.vars
    mgrid_dims = mgridnetcdf.dim

    nr = 0
    nz = 0
    nθ = 0
    external_coils = 1
    for i in 1:length(mgrid_dims.vals)
        if isassigned(mgrid_dims.vals, i)
            val = mgrid_dims.vals[i]
            if val.name == "rad"
                nr = Int(val.dimlen)
            elseif val.name == "zee"
                nz = Int(val.dimlen)
            elseif val.name == "phi"
                nθ = Int(val.dimlen)
            elseif val.name == "external_coils"
                external_coils = Int(val.dimlen)
            end
        end
    end
    rmin = NetCDF.readvar(mgrid_vars["rmin"])[1]
    rmax = NetCDF.readvar(mgrid_vars["rmax"])[1]
    zmin = NetCDF.readvar(mgrid_vars["zmin"])[1]
    zmax = NetCDF.readvar(mgrid_vars["zmax"])[1]
    nfp = NetCDF.readvar(mgrid_vars["nfp"])[1]
    
    r = range(rmin, rmax, nr)
    z = range(zmin, zmax, nz)
    #Δθ = (2*π/nfp)/(nθ + 1)
    θ = range(0, 2*π/nfp, nθ+1)

    fullSize = (length(r), length(θ), length(z))
    r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
    θ_grid = reshape(repeat(θ,inner=length(r),outer=length(z)),fullSize)
    z_grid = reshape(repeat(z,inner=length(r)*length(θ)),fullSize)
    field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))

    Br = permutedims(NetCDF.readvar(mgrid_vars["br_001"]).*current[1],[1,3,2])
    Bz = permutedims(NetCDF.readvar(mgrid_vars["bz_001"]).*current[1],[1,3,2])
    Bθ = permutedims(NetCDF.readvar(mgrid_vars["bp_001"]).*current[1],[1,3,2])

    for i in 2:external_coils
        if i > length(current)
            continue
        end
        numeral = lpad(i,3,'0')
        Br += permutedims(NetCDF.readvar(mgrid_vars["br_"*numeral]).*current[i],[1,3,2])
        Bz += permutedims(NetCDF.readvar(mgrid_vars["bz_"*numeral]).*current[i],[1,3,2])
        Bθ += permutedims(NetCDF.readvar(mgrid_vars["bp_"*numeral]).*current[i],[1,3,2])
    end

    Brf = zeros(nr, nθ+1, nz)
    Bzf = zeros(nr, nθ+1, nz)
    Bθf = zeros(nr, nθ+1, nz)

    Brf[:,1:end-1,:] = Br[:,:,:]
    Bzf[:,1:end-1,:] = Bz[:,:,:]
    Bθf[:,1:end-1,:] = Bθ[:,:,:]
    
    Brf[:,end,:] = Br[:,1,:]
    Bzf[:,end,:] = Bz[:,1,:]
    Bθf[:,end,:] = Bθ[:,1,:]

    #No potential, so return here
    if !("ar_001" in keys(mgridnetcdf.vars))
      return MagneticField(field_coords, Brf, Bθf, Bzf, nfp=nfp)
    end

    #Some mgrids will also have magnetic potentials, so let's add them also
    Ar = permutedims(NetCDF.readvar(mgrid_vars["ar_001"]).*current[1],[1,3,2])
    Az = permutedims(NetCDF.readvar(mgrid_vars["az_001"]).*current[1],[1,3,2])
    Aθ = permutedims(NetCDF.readvar(mgrid_vars["ap_001"]).*current[1],[1,3,2])

    for i in 2:external_coils
        if i > length(current)
            continue
        end
        numeral = lpad(i,3,'0')
        Ar += permutedims(NetCDF.readvar(mgrid_vars["ar_"*numeral]).*current[i],[1,3,2])
        Az += permutedims(NetCDF.readvar(mgrid_vars["az_"*numeral]).*current[i],[1,3,2])
        Aθ += permutedims(NetCDF.readvar(mgrid_vars["ap_"*numeral]).*current[i],[1,3,2])
    end

    Arf = zeros(nr, nθ+1, nz)
    Azf = zeros(nr, nθ+1, nz)
    Aθf = zeros(nr, nθ+1, nz)

    Arf[:,1:end-1,:] = Ar[:,:,:]
    Azf[:,1:end-1,:] = Az[:,:,:]
    Aθf[:,1:end-1,:] = Aθ[:,:,:]
    
    Arf[:,end,:] = Ar[:,1,:]
    Azf[:,end,:] = Az[:,1,:]
    Aθf[:,end,:] = Aθ[:,1,:]
    

    return MagneticField(field_coords, Brf, Bθf, Bzf, Arf, Aθf, Azf, nfp=nfp)

end
"""
  read_mgrid_h5(filename; currents=nothing)

This reads an mgrid file of the custom h5 type.  
If the currents keyword is set it will return a single
magnetic grid object.  Otherwise it will return a multiple grid object
"""
function read_mgrid_h5(filename::AbstractString; currents=nothing)
  file_id = h5open(filename)
  #count the number of groups
  ngroup = 0
  for group_id in file_id
    ngroup += 1
  end
  #there should always be at least two groups
  if ngroup < 2
    println("Not enough groups, bad h5 file")
    return nothing
  end

  # can make this fill in 0s eventually
  if currents != nothing && length(currents) != ngroup - 1
    println("Not enough/too many currents, need ",ngroup-1, "elements in currents file")
    return nothing
  end

  #get the governing data axes
  nfp = read(file_id["axes"]["nfp"])
  r = read(file_id["axes"]["r"])
  θ = read(file_id["axes"]["theta"])
  z = read(file_id["axes"]["z"])
  T = typeof(file_id["coil1"]["Br"][1,1,1])
  rrange = range(r[1],r[end],length(r))
  zrange = range(z[1],z[end],length(z))
  θrange = range(θ[1],θ[end],length(θ))

  #coords for making the grid (we repeat this code several times, todo: refactor it)
  fullSize = (length(rrange), length(θrange), length(zrange))
  r_grid = reshape(repeat(rrange,outer=length(zrange)*length(θrange)),fullSize)
  θ_grid = reshape(repeat(θrange,inner=length(rrange),outer=length(zrange)),fullSize)
  z_grid = reshape(repeat(zrange,inner=length(rrange)*length(θrange)),fullSize)
  field_coords = StructArray{Cylindrical}((r_grid, θ_grid, z_grid))
  if currents == nothing
    mgrid = MultipleFieldGrid(rrange, θrange, zrange, nfp, ngroup - 1)

  end


  Br = Array{T}(undef, length(r), length(θ), length(z))
  Bz = similar(Br)
  Bθ = similar(Br)
  Ar = similar(Br)
  Az = similar(Br)
  Aθ = similar(Br)
  Brt = similar(Br)
  Bzt = similar(Br)
  Bθt = similar(Br)
  Art = similar(Br)
  Azt = similar(Br)
  Aθt = similar(Br)

  for i in 1:ngroup-1
    group_name = "coil"*string(i)
    Brt .= read(file_id[group_name]["Br"])
    Bzt .= read(file_id[group_name]["Bz"])
    Bθt .= read(file_id[group_name]["Btheta"])
    Art .= read(file_id[group_name]["Ar"])
    Azt .= read(file_id[group_name]["Az"])
    Aθt .= read(file_id[group_name]["Atheta"])
    #
    if currents != nothing
      Br .= currents[i] .* Brt
      Bz .= currents[i] .* Bzt
      Bθ .= currents[i] .* Bθt
      Ar .= currents[i] .* Art
      Az .= currents[i] .* Azt
      Aθ .= currents[i] .* Aθt
    else
      mf = MagneticField(field_coords, Brt, Bθt, Bzt, Art, Aθt, Azt)
      mgrid.magnetic_field[i] = mf
    end
  end
  if currents != nothing
    return MagneticField(field_coords, Br, Bθ, Bz, Ar, Aθ, Az)
  else
    return mgrid
  end
  
end  
