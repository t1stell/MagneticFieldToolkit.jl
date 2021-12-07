
struct BFieldInterpolator{T<:AbstractFloat}
  Br_coeffs::Vector{Interpolations.Extrapolation}
  Bz_coeffs::Vector{Interpolations.Extrapolation}
  Bϕ_coeffs::Vector{Interpolations.Extrapolation}
  space::Chebyshev{ClosedInterval{T}, T}
end


mutable struct BField{T<:AbstractFloat}
  nr::Integer
  nz::Integer
  nphi::Integer
  nfp::Integer
  n_modes::Integer

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
                T::Type=Float64,
               )

  nfp = 1

  rmin = zero(T)
  zmin = zero(T)
  rmax = zero(T)
  zmax = zero(T)

  r = zeros(T,nr)
  z = zeros(T,nz)
  phi = zeros(T,nz)

  r_range = 0.0:1.0:1.0
  z_range = 0.0:1.0:1.0
  ϕ_range = 0.0:1.0:1.0

  #note this gets read in reverse order from the netcdf file
  br_grid = zeros(T,nr,nz,nphi)
  bz_grid = zeros(T,nr,nz,nphi)
  bp_grid = zeros(T,nr,nz,nphi)

  BField{T}(nr,nz,nphi,n_modes,nfp,rmin,zmin,rmax,zmax,
            r,z,phi,r_range, z_range, ϕ_range,
            br_grid, bz_grid, bp_grid)

end
