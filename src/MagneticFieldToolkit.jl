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
using PlasmaEquilibriumToolkit


# File load 
export read_mgrid, read_bmw, read_vmec_coils, read_mgrid_h5

#Structs and constructors
export CoilSet

# Field calculations
export compute_magnetic_potential, compute_magnetic_field
export generate_mgrid

# Following calculations
export follow_field, follow_field_s

include("bfieldUtils.jl")
include("ReadMagneticField.jl")
include("QuadraticFluxMinimize.jl")
include("Coils.jl")
include("followField.jl")

function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("plot_makie.jl")
end

end
