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
export read_mgrid, read_bmw, read_vmec_coils

# Field calculations
export compute_magnetic_potential, compute_magnetic_field
export generate_mgrid

include("bfieldUtils.jl")
include("ReadMagneticField.jl")
include("followField.jl")
include("QuadraticFluxMinimize.jl")
include("Coils.jl")

function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("plot_makie.jl")
end

end
