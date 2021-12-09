module FieldlineTracing

using StaticArrays
using NetCDF
using Interpolations
using Polyester
using DifferentialEquations

include("bfieldUtils.jl")
include("bfieldTypes.jl")
include("bfield.jl")
include("followField.jl")

end
