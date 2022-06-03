module MagneticFieldToolkit
using Requires
using StaticArrays
using NetCDF
using Interpolations
using Polyester
using DifferentialEquations

include("bfieldUtils.jl")
include("bfieldTypes.jl")
include("bfield.jl")
include("followField.jl")

function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("plot_makie.jl")
end

end
