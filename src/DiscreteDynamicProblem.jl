

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
struct DiscreteDynamicProblem{nStateVars,nChoiceVars,C,G,IP,IF,nChoiceVarsZero, C0} <: DDM

    β::Real # important to specify here for type stability

    # I think it should be possible to specify the number of inputs of the function
    rewardfunc::Function
    transfunc::Function
    intdim::Type{ID} where ID<:DDMIntDim # not really necessary, can just have it as a parameter in type...

    tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{nChoiceVars, C}

    # distribution of shocks
    # allows few univariate only normal multivariate for now
    shockdist::Distribution  # for smm it's important though that the shock distribution stays the same, does not
    # depend on parameters!!!

    # which state variables are endogenous
    # bEndogStateVars::E

    """The partial reward function."""
    grossprofits::G

    """Objective function at t=0 to choose initial endogenous state variables."""
    initializationproblem::IP

    """Function that maps initial shocks and initial policy into state variables."""
    initializefunc::IF

    tChoiceVectorsZero::NTuple{nChoiceVarsZero, C0}

end


function createDiscreteDynamicProblem(
            # params::ModelParams,
            β::Real,
            rewardfunc::Function,
            transfunc::Function,
            tStateVectors::NTuple{S, Vector{Float64}},
            tChoiceVectors::NTuple{C, typeC},
            shockdist::Distribution;
            # vWeights::Vector{Float64},
            # mShocks::Array{Float64,2};
            intdim::Symbol = :All,
            # bEndogStateVars::Union{Vector{Bool},Nothing} = nothing,
            grossprofits::Union{Function,Nothing} = nothing,
            initializationproblem::Union{Function,Nothing} = nothing,
            initializefunc::Union{Function,Nothing} = nothing,
            tChoiceVectorsZero::NTuple{C0, typeC0} = tChoiceVectors,
            ) where {I <: DDMIntDim, S, C, typeC<:Union{Vector{Float64}, Int64},
                C0, typeC0<:Union{Vector{Float64}, Int64}}

    # can do stuff that the user does not interact with
    # nStateVars = length(tStateVectors)

    # tChoiceVectors =
    #
    # if typeof(tChoiceVectors) <: Vector{Bool}
    #     nChoiceVars = sum(tChoiceVectors)
    # else
    #     nChoiceVars = length(tStateVectors)
    # end

    DiscreteDynamicProblem(
        β,
        rewardfunc,
        transfunc,
        eval(intdim),
        tStateVectors,
        tChoiceVectors,
        shockdist,
        # vWeights,
        # mShocks,
        # bEndogStateVars,
        grossprofits,
        initializationproblem,
        initializefunc,
        tChoiceVectorsZero,
        )
end

const DDP = DiscreteDynamicProblem

separable(p::DDP) = typeof(p.tChoiceVectors) <: NTuple{N, Int} where N
