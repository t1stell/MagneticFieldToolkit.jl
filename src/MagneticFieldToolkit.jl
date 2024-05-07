module MagneticFieldToolkit
using Requires
using StaticArrays
using StructArrays
using NetCDF
using Interpolations
using Polyester
using OrdinaryDiffEq
using LinearAlgebra
using HCubature
using CoordinateTransformations
using PolygonOps
using PlasmaEquilibriumToolkit
using StellaratorGrids


# File load 
export read_mgrid, read_bmw, read_vmec_coils, read_mgrid_h5, read_coil_files

#Structs and constructors
export CoilSet

# Field calculations
export compute_magnetic_potential, compute_magnetic_field
export generate_mgrid, extreme_coils

# Following calculations
export follow_field, follow_field_s, follow_to_wall

#Toroidal Flux
export toroidal_flux, iota_from_poincare, iota_at_r

include("ReadMagneticField.jl")
include("QuadraticFluxMinimize.jl")
include("Coils.jl")
include("bfieldUtils.jl")
include("followField.jl")
include("ToroidalFlux.jl")

function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("plot_makie.jl")
end

end
