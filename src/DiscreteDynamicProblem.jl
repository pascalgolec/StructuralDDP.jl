

abstract type ModelParams end

struct DiscreteDynamicProblem <: DDM # DiscreteDynamicModel

    rewardfunc::Function
    transfunc::Function
    grossprofits::Function

    initializationproblem::Union{Function, Nothing}
    initializefunc::Union{Function, Nothing}

    # params::NeoClassicalSimpleParams{Float64, Int64} # important to specify here for type stability
    params::ModelParams # important to specify here for type stability

    tStateVectors::NTuple{N, Vector{Float64}} where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::Union{Tuple{Vector{Float64}},Tuple{Vector{Float64},Vector{Float64}}}

    # which state variables are endogenous
    bEndogStateVars::Vector{Bool}

    # quadrature for calculating expectations
    vWeights::Vector{Float64}
    mShocks::Array{Float64,2}

end
