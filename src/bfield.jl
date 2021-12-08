
#read an mgrid file and create a relevant struct

#function readMgrid(filename::String)
#  mgrid = NetCDF.open(mgridfile)
#

#read a bmw file and create a relevant field struct
"""
    read_BMW(filename; n_modes=0, padding=5)

Read a NetCDF BMW file and store the result in a `BField` type.

# Arguemnts:
- `filename::AbstractString`   # NetCDF BMW file name
- `n_modes::Integer` = 0   # Number of modes used for Chebyshev interpolation
- `padding::Integer` = 5   # Number of slices in ϕ padded to both ends of the magnetic grid arrays for accurate interpolation

# See also: [`BField`](@ref)
"""
function readBMW(filename::AbstractString;
                  n_modes::Integer=0,
                  padding::Integer=5,
                 )
  bmwnetcdf = NetCDF.open(filename)
  bmwVars = bmwnetcdf.vars
  bmwDims = bmwnetcdf.dim
#
#
  #there's probably a better way to do this, but damn the docs are awful
  nr = 0
  nz = 0
  nphi = 0
  for i in 1:length(bmwDims.vals)
    #why the fuck there are undefined values here, who knows???
    if isassigned(bmwDims.vals, i)
      val = bmwDims.vals[i]
      if val.name == "r"
        nr = Int(val.dimlen)
        @debug "Size of grid in r direction: $nr"
      elseif val.name == "z"
        nz = Int(val.dimlen)
        @debug "Size of grid in z direction: $nz"
      elseif val.name == "phi"
        nphi = Int(val.dimlen)
        @debug "Size of grid in ϕ direction: $nphi"
      end
    end
  end

  #How many values are in the FFT
  n_modes = n_modes > div(nphi,2) ? n_modes : div(nphi,2)
  @debug "Number of modes used for Chebyshev interpolation: $n_modes"
  floatType = eltype(NetCDF.readvar(bmwVars["rmin"]))

  bmw = BField(nr, nz, nphi; n_modes=n_modes, padding=padding, T=floatType)
  bmw.nr = nr
  bmw.nz = nz
  bmw.nphi = nphi
  bmw.n_modes = n_modes
  bmw.padding = padding

  extended_data = zeros(floatType, nr, nz, nphi+2*padding)
  #this loop looks through all the variables in bmwVars and
  #checks if any appear in the array, if they do, then it
  #copies the data into them
  for var in bmwVars
    varSymbol = Symbol(var.first)
    if hasproperty(bmw, varSymbol)
      T = typeof(getfield(bmw, varSymbol))
      var_data = NetCDF.readvar(var.second)
      if T <: Array
        # Pad either end of the array using periodic BCs
        extended_data[:, :, 1+padding:nphi+padding] = convert.(eltype(T), var_data)
        extended_data[:, :, 1:padding] = extended_data[:, :, nphi+1:nphi+padding]
        extended_data[:, :, nphi+padding+1:end] = extended_data[:, :, padding+1:2*padding]
        setfield!(bmw, varSymbol, copy(extended_data))
        @debug "Set BField variable $varSymbol"
        fill!(extended_data, zero(eltype(T)))
      elseif T <: Number
        setfield!(bmw, varSymbol,convert(T,var_data[1]))
        @debug "Set BField variable $varSymbol = $(getfield(bmw, varSymbol))"
      elseif T <: String
        setfield!(bmw, varSymbol,string(strip(string(var_data...))))
        @debug "Set BField variable $varSymbol = $(getfield(bmw, varSymbol))"
      else
        throw(TypeError("No type given for $(var)!"))
      end
    end
  end

  phimax = (2π)/bmw.nfp
  bmw.r = LinRange(bmw.rmin, bmw.rmax, bmw.nr)
  bmw.z = LinRange(bmw.zmin, bmw.zmax, bmw.nz)
  bmw.phi = LinRange(0, phimax, bmw.nphi)

  bmw.r_range = vector2range(bmw.r)
  bmw.z_range = vector2range(bmw.z)
  bmw.ϕ_range = vector2range(bmw.phi)

  return bmw
end

function BFieldInterpolator(bfield::BField{T};
                            coeff_threshold::T = 10^(-4),
                           ) where {T}
  ϕ_knots = range(first(bfield.ϕ_range)-step(bfield.ϕ_range)*bfield.padding,step=step(bfield.ϕ_range),length=length(bfield.ϕ_range)+2*bfield.padding)
  r_knots = bfield.r_range
  z_knots = bfield.z_range
  nr = bfield.nr
  nz = bfield.nz

  N = length(ϕ_knots)
  M = bfield.n_modes

  S = Chebyshev(first(ϕ_knots)..last(ϕ_knots))
  @debug "Interpolation space: $S"
  V = Array{T}(undef, N, M)
  @debug "Size of Vandermonde matrix: $(size(V))"

  @batch minbatch=64 for k in 1:M
    V[:, k] = Fun(S, [zeros(k-1);1]).(ϕ_knots)
  end
  @debug "Constructed Vandermonde matrix"

  Br_coeff_grid = Array{Float64,3}(undef, M, nr, nz)
  Bz_coeff_grid = similar(Br_coeff_grid)
  Bϕ_coeff_grid = similar(Br_coeff_grid)
  
  for z in 1:nz
    for r in 1:nr
      @inbounds Br_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.br_grid, r, z, :)))
      @inbounds Bz_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.bz_grid, r, z, :)))
      @inbounds Bϕ_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.bp_grid, r, z, :)))
    end
  end
  @debug "Constructed coefficient grids"


  Br_coeffs = Vector{Interpolations.Extrapolation}(undef,M)
  Bz_coeffs = similar(Br_coeffs)
  Bϕ_coeffs = similar(Br_coeffs)
  @batch for m in 1:M
    Br_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Br_coeff_grid[m,:,:])
    Bz_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Bz_coeff_grid[m,:,:])
    Bϕ_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Bϕ_coeff_grid[m,:,:])
  end
  @debug "Constructed coefficient interpolators"

  BFieldInterpolator{T}(Br_coeffs, Bz_coeffs, Bϕ_coeffs, S, 1:bfield.n_modes, bfield.nfp)

end

function readMgrid(filename::AbstractString, current::Vector{Float64};
                   n_modes::Integer = 0, padding::Integer=5,
                  )
  
  mgridnetcdf = NetCDF.open(filename)
  mgridVars = mgridnetcdf.vars
  mgridDims = mgridnetcdf.dim

  nr = 0
  nz = 0
  nphi = 0
  external_coils = 1
  for i in 1:length(mgridDims.vals)
    if isassigned(mgridDims.vals, i)
      val = mgridDims.vals[i]
      if val.name == "rad"
        nr = Int(val.dimlen)
      elseif val.name == "zee"
        nz = Int(val.dimlen)
      elseif val.name == "phi"
        nphi = Int(val.dimlen)
      elseif val.name == "external_coils"
        external_coils == Int(val.dimlen)
      end
    end
  end
  n_modes = n_modes > div(nphi, 2) ? n_modes : div(nphi, 2)
  floatType = eltype(NetCDF.readvar(mgridVars["rmin"]))

  mg = BField(nr, nz, nphi; n_modes=n_modes, 
               external_coils = external_coils, padding=padding) 

  mg.nr = nr
  mg.nz = nz
  mg.nphi = nphi
  mg.n_modes = n_modes
  mg.padding = padding
  mg.external_coils = external_coils

  bp_grid = zeros(floatType, nr, nz, nphi + 2*padding)
  br_grid = zeros(floatType, nr, nz, nphi + 2*padding)
  bz_grid = zeros(floatType, nr, nz, nphi + 2*padding)

  b_temp = zeros(floatType, nr, nz, nphi + 2*padding)

  for var in mgridVars
    varSymbol = Symbol(var.first)
    var_data = NetCDF.readvar(var.second)
    varString = String(var.first)
   
    #Special handling for the field values
    if  (occursin("bp_", varString) || occursin("bz_", varString) || 
        occursin("br_", varString))

      
      T = typeof(getfield(mg, Symbol("bp_grid")))
      #get the index
      i = parse(Int64, varString[4:end])
      ci = current[i]
      fill!(b_temp, zero(eltype(T)))

      b_temp[:,:,1+padding:nphi+padding] = convert.(eltype(T), var_data)
      b_temp[:,:,1:padding] = b_temp[:, :, nphi+1:nphi+padding]
      b_temp[:,:,nphi+padding+1:end]=b_temp[:,:,padding+1:2*padding]
      b_temp *= ci
      if occursin("bp_", varString)
        bp_grid += copy(b_temp)
      elseif occursin("bz_", varString)
        bz_grid += copy(b_temp)
      elseif occursin("br_", varString)
        br_grid += copy(b_temp)
      end
    
    #handle the other properties which are named the same
    elseif hasproperty(mg, varSymbol)
      T = typeof(getfield(mg, varSymbol))
      if T <: Number
        setfield!(mg, varSymbol, convert(T, var_data[1]))
      elseif T <: String
        setfield!(mg, varSymbol, string(strip(string(var_data...))))
      else
        throw(TypeError("No type given for $(var)!"))
      end
    end

  end

  setfield!(mg, Symbol("bp_grid"), copy(bp_grid))
  setfield!(mg, Symbol("br_grid"), copy(br_grid))
  setfield!(mg, Symbol("bz_grid"), copy(bz_grid))


  phimax = (2π)/mg.nfp
  mg.r = LinRange(mg.rmin, mg.rmax, mg.nr)
  mg.z = LinRange(mg.zmin, mg.zmax, mg.nz)
  mg.phi = LinRange(0, phimax, mg.nphi)

  mg.r_range = vector2range(mg.r)
  mg.z_range = vector2range(mg.z)
  mg.ϕ_range = vector2range(mg.phi)

  return mg

end  


function (itp::BFieldInterpolator)(r::T,
                                   z::T,
                                   ϕ::T;
                                   xyz = false
                                  ) where {T}
  ϕ = mod(ϕ, 2*π/itp.nfp)
  Br = Fun(itp.space, [itp.Br_coeffs[m](r,z) for m in itp.modes])
  Bz = Fun(itp.space, [itp.Bz_coeffs[m](r,z) for m in itp.modes])
  Bϕ = Fun(itp.space, [itp.Bϕ_coeffs[m](r,z) for m in itp.modes])

  if xyz
    Br = Br(ϕ)
    Bϕ = Bϕ(ϕ)
    Bx = Br * cos(ϕ) - Bϕ * sin(ϕ)
    By = Br * sin(ϕ) + Bϕ * cos(ϕ)
    return (Bx, By, Bz(ϕ))
  end
  return (Br(ϕ), Bz(ϕ), Bϕ(ϕ))
end
