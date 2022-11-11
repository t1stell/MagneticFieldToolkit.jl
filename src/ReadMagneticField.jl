
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

