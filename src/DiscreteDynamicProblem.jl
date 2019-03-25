

abstract type ModelParams end

struct DiscreteDynamicProblem <: SingleChoiceVar # DiscreteDynamicModel

    rewardfunc::Function
    transfunc::Function
    grossprofits::Union{Function, Nothing}

    # params::NeoClassicalSimpleParams{Float64, Int64} # important to specify here for type stability
    params::ModelParams # important to specify here for type stability

    tStateVectors::NTuple{2, Vector{Float64}} # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{1, Vector{Float64}}

    # which state variables are endogenous
    bEndogStateVars::Vector{Bool}

    # quadrature for calculating expectations
    vWeights::Vector{Float64}
    mShocks::Array{Float64,2}

end
