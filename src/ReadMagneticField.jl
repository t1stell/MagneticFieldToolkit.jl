
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
    θ = range(0, 2*π/nfp, nθ)

    fullSize = (length(r), length(z), length(θ))
    r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
    z_grid = reshape(repeat(z,inner=length(r),outer=length(θ)),fullSize)
    θ_grid = reshape(repeat(θ,inner=length(z)*length(r)),fullSize)
    field_coords = StructArray{Cylindrical}((r_grid, z_grid, θ_grid))

    Br = NetCDF.readvar(bmw_vars["br_grid"])
    Bz = NetCDF.readvar(bmw_vars["bz_grid"])
    Bθ = NetCDF.readvar(bmw_vars["bp_grid"])

    Ar = NetCDF.readvar(bmw_vars["ar_grid"])
    Az = NetCDF.readvar(bmw_vars["az_grid"])
    Aθ = NetCDF.readvar(bmw_vars["ap_grid"])
    
    return MagneticField(field_coords, Br, Bz, Bθ, Ar, Az, Aθ, nfp=nfp)

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
    θ = range(0, 2*π/nfp, nθ)

    fullSize = (length(r), length(z), length(θ))
    r_grid = reshape(repeat(r,outer=length(z)*length(θ)),fullSize)
    z_grid = reshape(repeat(z,inner=length(r),outer=length(θ)),fullSize)
    θ_grid = reshape(repeat(θ,inner=length(z)*length(r)),fullSize)
    field_coords = StructArray{Cylindrical}((r_grid, z_grid, θ_grid))

    Br = NetCDF.readvar(mgrid_vars["br_001"]).*current[1]
    Bz = NetCDF.readvar(mgrid_vars["bz_001"]).*current[1]
    Bθ = NetCDF.readvar(mgrid_vars["bp_001"]).*current[1]

    for i in 2:external_coils
        if i > length(current)
            continue
        end
        numeral = lpad(i,3,'0')
        Br += NetCDF.readvar(mgrid_vars["br_"*numeral]).*current[i]
        Bz += NetCDF.readvar(mgrid_vars["bz_"*numeral]).*current[i]
        Bθ += NetCDF.readvar(mgrid_vars["bp_"*numeral]).*current[i]
    end


    return MagneticField(field_coords, Br, Bz, Bθ, nfp=nfp)

end

