using DiscreteDynamicProgramming
using Test
mytol = 1e-4

# # load models
# include("models/CapitalAdjustModel.jl")
# include("models/CapitalAdjustModel2.jl")
# createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

include("test_transitionmatrix.jl")
include("test_rewardmatrix_options.jl")
include("test_intdim.jl")
include("test_monotonicity_concavity.jl")
include("test_simulate.jl")
# include("test_speed.jl")
