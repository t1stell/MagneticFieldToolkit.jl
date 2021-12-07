module FieldlineTracing

using StaticArrays
using NetCDF
using Interpolations
using ApproxFun
using IntervalSets
using DifferentialEquations
using Polyester

include("bfieldUtils.jl")
include("bfieldTypes.jl")
include("bfield.jl")
include("followField.jl")

end
