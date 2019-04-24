using DiscreteDynamicProgramming
using Test
mytol = 1e-4

# load models
include("models/Neoclassical_user.jl")
include("models/Intangible_user.jl")
createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

include("test_transitionmatrix.jl")
include("test_rewardmatrix_options.jl")
include("test_intdim.jl")
include("test_monotonicity_concavity.jl")
include("test_simulate.jl")
include("test_speed.jl")
