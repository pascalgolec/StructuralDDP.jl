

abstract type ModelParams end


struct DiscreteDynamicProblem{nStateVars, nChoiceVars, Intdim} <: DDM where K<:DDMIntDim # DiscreteDynamicModel
    # with parameters, we create a family of types
    # can also add where F <: Union{Function, Nothing} later

    # I think it should be possible to specify the number of inputs of the function
    rewardfunc::Function
    transfunc::Function
    grossprofits::Function

    # initializationproblem::Union{Function, Nothing} # this is bad, because type will always stay ambiguous!!
    # need to parametrize above
    # initializefunc::Union{Function, Nothing}
    initializationproblem::Function
    initializefunc::Function

    # params::NeoClassicalSimpleParams{Float64, Int64} # important to specify here for type stability
    params::ModelParams # important to specify here for type stability

    intdim::Type{Intdim} # not really necessary, can just have it as a parameter in type...

    # tStateVectors::NTuple{N, Vector{Float64}} where N # can use NTuple{N, Vector{Float64}} where N
    # tChoiceVectors::Union{Tuple{Vector{Float64}},Tuple{Vector{Float64},Vector{Float64}}}

    tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{nChoiceVars, Vector{Float64}}

    # which state variables are endogenous
    bEndogStateVars::Vector{Bool}

    # quadrature for calculating expectations
    vWeights::Vector{Float64}
    mShocks::Array{Float64,2}

end
