using PlasmaEquilibriumToolkit
using StaticArrays
using StructArrays
using CoordinateTransformations
using StellaratorGrids
using MagneticFieldToolkit
using Test

@testset "MagneticFieldToolkit.jl" begin
    include("file_load_tests.jl")
    include("follow_field_test.jl")
    include("coil_tests.jl")
end
