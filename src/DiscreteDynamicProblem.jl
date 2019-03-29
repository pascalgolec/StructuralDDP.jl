

abstract type ModelParams end

#
# struct DiscreteDynamicProblem{nStateVars, nChoiceVars, Intdim} <: DDM where K<:DDMIntDim # DiscreteDynamicModel
#     # with parameters, we create a family of types
#     # can also add where F <: Union{Function, Nothing} later
#
#     # I think it should be possible to specify the number of inputs of the function
#     rewardfunc::Function
#     transfunc::Function
#
#     """The partial reward function."""
#     grossprofits::Function
#
#     # initializationproblem::Union{Function, Nothing} # this is bad, because type will always stay ambiguous!!
#     # need to parametrize above
#     # initializefunc::Union{Function, Nothing}
#     initializationproblem::Function
#     initializefunc::Function
#
#     # params::NeoClassicalSimpleParams{Float64, Int64} # important to specify here for type stability
#     params::ModelParams # important to specify here for type stability
#
#     intdim::Type{Intdim} # not really necessary, can just have it as a parameter in type...
#
#     # tStateVectors::NTuple{N, Vector{Float64}} where N # can use NTuple{N, Vector{Float64}} where N
#     # tChoiceVectors::Union{Tuple{Vector{Float64}},Tuple{Vector{Float64},Vector{Float64}}}
#
#     tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
#     tChoiceVectors::NTuple{nChoiceVars, Vector{Float64}}
#
#     # which state variables are endogenous
#     bEndogStateVars::Vector{Bool}
#
#     # quadrature for calculating expectations
#     vWeights::Vector{Float64}
#     mShocks::Array{Float64,2}
#
# end


# <: DDM
struct DiscreteDynamicProblem{nStateVars,nChoiceVars,E,G,IP,IF} <: DDM
# struct DiscreteDynamicProblem#{E,G,IP,IF}
    # with parameters, we create a family of types
    # can also add where F <: Union{Function, Nothing} later

    # params::NeoClassicalSimpleParams{Float64, Int64} # important to specify here for type stability
    """The parameters of the model are inside the structure params. (which could change in the future)."""
    params::ModelParams # important to specify here for type stability

    # I think it should be possible to specify the number of inputs of the function
    rewardfunc::Function
    transfunc::Function

    intdim::Type{ID} where ID<:DDMIntDim # not really necessary, can just have it as a parameter in type...

    tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{nChoiceVars, Vector{Float64}}

    # distribution of shocks
    # transitionmatrix() allows few univariate only normal multivariate for now
    shockdist::Distribution  # for smm it's important though that the shock distribution stays the same, does not
    # depend on parameters!!!

    # which state variables are endogenous
    bEndogStateVars::E

    """The partial reward function."""
    grossprofits::G

    # initializationproblem::Union{Function, Nothing} # this is bad, because type will always stay ambiguous!!
    # need to parametrize above
    initializationproblem::IP
    initializefunc::IF

end


function DiscreteDynamicProblem(
            params::ModelParams,
            rewardfunc::Function,
            transfunc::Function,
            intdim::Type{I},
            tStateVectors::NTuple{S, Vector{Float64}},
            tChoiceVectors::NTuple{C, Vector{Float64}},
            shockdist::Distribution;
            # vWeights::Vector{Float64},
            # mShocks::Array{Float64,2};
            bEndogStateVars::Union{Vector{Bool},Nothing} = nothing,
            grossprofits::Union{Function,Nothing} = nothing,
            initializationproblem::Union{Function,Nothing} = nothing,
            initializefunc::Union{Function,Nothing} = nothing,
            ) where {I <: DDMIntDim, S, C}

    # can do stuff that the user does not interact with

    DiscreteDynamicProblem(params,
        rewardfunc,
        transfunc,
        intdim,
        tStateVectors,
        tChoiceVectors,
        shockdist,
        # vWeights,
        # mShocks,
        bEndogStateVars,
        grossprofits,
        initializationproblem,
        initializefunc,
        )
end
