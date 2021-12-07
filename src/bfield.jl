
#read an mgrid file and create a relevant struct

#function readMgrid(filename::String)
#  mgrid = NetCDF.open(mgridfile)
#

#read a bmw file and create a relevant field struct
function read_BMW(filename::AbstractString;
                  n_modes::Integer=0,
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

  bmw = BField(nr, nz, nphi; n_modes=n_modes, T=floatType)
  bmw.nr = nr
  bmw.nz = nz
  bmw.nphi = nphi
  bmw.n_modes = n_modes

  #this loop looks through all the variables in bmwVars and
  #checks if any appear in the array, if they do, then it
  #copies the data into them
  for var in bmwVars
    varSymbol = Symbol(var.first)
    if hasproperty(bmw, varSymbol)
      T = typeof(getfield(bmw, varSymbol))
      var_data = NetCDF.readvar(var.second)
      if T <: Array
        setfield!(bmw, varSymbol, convert.(eltype(T), var_data))
        @debug "Set BField variable $varSymbol"
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
  ϕ_knots = bfield.ϕ_range
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
    @batch for r in 1:nr
      Br_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.br_grid, r, z, :)))
      Bz_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.bz_grid, r, z, :)))
      Bϕ_coeff_grid[:, r, z] = coefficients(Fun(S, V \ view(bfield.bp_grid, r, z, :)))
    end
  end


  Br_coeffs = Vector{Interpolations.Extrapolation}(undef,M)
  Bz_coeffs = similar(Br_coeffs)
  Bϕ_coeffs = similar(Br_coeffs)
  @batch for m in 1:M
    Br_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Br_coeff_grid[m,:,:])
    Bz_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Bz_coeff_grid[m,:,:])
    Bϕ_coeffs[m] = CubicSplineInterpolation((r_knots, z_knots), Bϕ_coeff_grid[m,:,:])
  end

  BFieldInterpolator{T}(Br_coeffs, Bz_coeffs, Bϕ_coeffs, S)

end

function (itp::BFieldInterpolator)(r::T,
                                   z::T,
                                   ϕ::T;
                                  ) where {T}
  modes = 1:length(itp.Br_coeffs)
  Br = Fun(itp.space, [itp.Br_coeffs[m](r,z) for m in modes])
  Bz = Fun(itp.space, [itp.Bz_coeffs[m](r,z) for m in modes])
  Bϕ = Fun(itp.space, [itp.Bϕ_coeffs[m](r,z) for m in modes])
  return (Br(ϕ), Bz(ϕ), Bϕ(ϕ))
end
