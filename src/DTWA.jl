module DTWA
export greet

using DifferentialEquations
using LinearAlgebra
using Optim
include("solve.jl")
include("magnetisation.jl")
include("spinsqueezingparam.jl")
include("3d.jl")

function greet(string)
    return string
end

end