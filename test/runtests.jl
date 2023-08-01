using MagneticFieldToolkit
using StaticArrays
using CoordinateTransformations
using StellaratorGrids
using Test

@testset "MagneticFieldToolkit.jl" begin
  include("follow_field_test.jl")
  include("coil_tests.jl")
end
