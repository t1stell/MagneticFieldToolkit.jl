
struct BFieldInterpolator{T<:AbstractFloat}
  Br_coeffs::Vector{Interpolations.Extrapolation}
  Bz_coeffs::Vector{Interpolations.Extrapolation}
  Bϕ_coeffs::Vector{Interpolations.Extrapolation}
  space::Space
  modes::UnitRange
end


mutable struct BField{T<:AbstractFloat}
  nr::Integer
  nz::Integer
  nphi::Integer
  nfp::Integer
  n_modes::Integer
  padding::Integer

  rmin::T
  zmin::T
  rmax::T
  zmax::T

  r::Vector{T}
  z::Vector{T}
  phi::Vector{T}
  r_range::StepRangeLen
  z_range::StepRangeLen
  ϕ_range::StepRangeLen

  br_grid::Array{T,3}
  bz_grid::Array{T,3}
  bp_grid::Array{T,3}
end

function BField(nr::Integer,
                nz::Integer,
                nphi::Integer;
                n_modes::Integer = div(nϕ,2),
                padding::Integer = 5,
                T::Type=Float64,
               )

  nfp = 1

  rmin = zero(T)
  zmin = zero(T)
  rmax = zero(T)
  zmax = zero(T)

  r = zeros(T,nr)
  z = zeros(T,nz)
  phi = zeros(T,nphi)

  r_range = 0.0:1.0:1.0
  z_range = 0.0:1.0:1.0
  ϕ_range = 0.0:1.0:1.0

  #note this gets read in reverse order from the netcdf file
  br_grid = zeros(T,nr,nz,nphi+2*padding)
  bz_grid = zeros(T,nr,nz,nphi+2*padding)
  bp_grid = zeros(T,nr,nz,nphi+2*padding)

  BField{T}(nr, nz, nphi, nfp, n_modes, padding,
            rmin, zmin, rmax, zmax,
            r, z, phi, r_range, z_range, ϕ_range,
            br_grid, bz_grid, bp_grid)

end
