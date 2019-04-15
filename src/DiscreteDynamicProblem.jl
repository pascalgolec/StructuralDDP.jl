
"""
$(TYPEDEF)

Defines an infinite-horizon discrete choice dynamic optimization problem.

# Fields

$(FIELDS)
"""
struct DiscreteDynamicProblem{nStateVars,nChoiceVars,C,G,IP,IF,nChoiceVarsZero, C0} <: DDM

	"""The discount factor for future rewards."""
    β::Float64 # important to specify here for type stability

    """The reward function defining current period rewards as a function of states and choices."""
    rewardfunc::Function

	"""The transition function of the state variables, depending on the
	integration dimension as a function of states, choices and shocks."""
    transfunc::Function

	"""The integration dimension of the transition function."""
    intdim::Type{ID} where ID<:DDMIntDim # not really necessary, can just have it as a parameter in type...

    tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{nChoiceVars, C}

    # distribution of shocks
    # allows few univariate only normal multivariate for now
    shockdist::Distribution  # for smm it's important though that the shock distribution stays the same, does not
    # depend on parameters!!!

    """The partial reward function."""
    grossprofits::G

    """Objective function at t=0 to choose initial endogenous state variables."""
    initializationproblem::IP

    """Function that maps initial shocks and initial policy into state variables."""
    initializefunc::IF

    tChoiceVectorsZero::NTuple{nChoiceVarsZero, C0}

end

const DDP = DiscreteDynamicProblem

function createDDP(
            # params::ModelParams,
            β::Real,
            rewardfunc::Function,
            transfunc::Function,
            tStateVectors::NTuple{dimStates, Vector{Float64}},
            tChoiceVectors::NTuple{dimChoices, typeC},
            shockdist::Distribution;
            # vWeights::Vector{Float64},
            # mShocks::Array{Float64,2};
            intdim::Symbol = :All,
            # bEndogStateVars::Union{Vector{Bool},Nothing} = nothing,
            grossprofits::Union{Function,Nothing} = nothing,
            initializationproblem::Union{Function,Nothing} = nothing,
            initializefunc::Union{Function,Nothing} = nothing,
            tChoiceVectorsZero::Union{NTuple{C0,Int64},Nothing} = nothing,
            ) where {I <: DDMIntDim, dimStates, dimChoices, typeC<:Union{Vector{Float64}, Int64},
                C0}

    # can do stuff that the user does not interact with
    # nStateVars = length(tStateVectors)

	intdim in (:All, :Separable, :Separable_ExogStates, :Separable_States) ||
		error("Provided wrong integration dimension $intdim.")

	!(eval(intdim) <: Separable_Union && typeC <: AbstractVector) || error(
	"Provide a tuple of integers pointing towards the state variables instead of redefining the choice variables.")

	# solver only supports correct order of endogenous state variables
	!(eval(intdim) <: Separable_Union && tChoiceVectors[1] != 1) || error(
	"Bad tChoiceVectors: the first state variable must be the (first) choice variable.")
	!(eval(intdim) <: Separable_Union && dimChoices==2 && tChoiceVectors[2] != 2) || error(
	"tChoiceVectors: the second state variable must be the second choice variable.")

	!(typeof(tChoiceVectorsZero) <: NTuple{C0,Int64} where C0 && tChoiceVectorsZero[1] != 1) || error(
	"Bad tChoiceVectorsZero: the first state variable must be the (first) choice variable in the intialization problem.")


    DiscreteDynamicProblem(
        β,
        rewardfunc,
        transfunc,
        eval(intdim),
        tStateVectors,
        tChoiceVectors,
        shockdist,
        grossprofits,
        initializationproblem,
        initializefunc,
        tChoiceVectorsZero,
        )
end



separable(p::DDP) = typeof(p.tChoiceVectors) <: NTuple{N, Int} where N

"""Retreive the state vectors that are not choice vectors from the state vector tuple."""
getnonchoicevars(p::DDM) = getnonchoicevars(p.tStateVectors, p.tChoiceVectors)
function getnonchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC}
	allvars = tuple(1:NS...)
	nonchoicevars = tuple(setdiff(Set(allvars), Set(tchoicevars))...)
	return getindex(tStateVectors, collect(nonchoicevars))
end
getnonchoicevars(p::DDP{dimS,dimC,C}) where {dimS,dimC,C<:AbstractVector} = tuple() # i.e. intdim = :All
getnonchoicevars(tStateVectors, tChoiceVectors::Nothing) = nothing
getnonchoicevarszero(p::DDM) = getnonchoicevars(p.tStateVectors, p.tChoiceVectorsZero)


"""Retreive the choice vectors from the state vector tuple."""
getchoicevars(p::DDM) = getchoicevars(p.tStateVectors, p.tChoiceVectors)
getchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC} =
	getindex(tStateVectors, collect(tchoicevars))
"""Retreive nothing if the choice vectors are provided."""
getchoicevars(tStateVectors::NTuple{N1,T1}, tChoiceVectors::NTuple{N2,Vector{T2}}) where
	{N1, T1, N2, T2} = tChoiceVectors
getchoicevars(tStateVectors, tChoiceVectors::Nothing) = nothing
"""Retreive the choice vectors to find policy at t=0."""
getchoicevarszero(p::DDM) = getchoicevars(p.tStateVectors, p.tChoiceVectorsZero)
