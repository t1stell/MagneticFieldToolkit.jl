using MagneticFieldToolkit
using StaticArrays
using CoordinateTransformations
using StellaratorGrids
using Test

@testset "MagneticFieldToolkit.jl" begin
    include("file_load_tests.jl")
    include("follow_field_test.jl")
    include("coil_tests.jl")
end
