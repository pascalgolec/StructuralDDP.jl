
# model definition and parameters
# abstract type DiscreteDynamicModel end
# const DDP = DiscreteDynamicModel

"""The integration dimension determines the in- and outputs of the transition function."""
abstract type IntegrationDimension end
const IntDim = IntegrationDimension
abstract type All <: IntDim end
abstract type Separable <: IntDim end
abstract type Separable_States <: IntDim end
abstract type Separable_ExogStates <: IntDim end
const Separable_Union = Union{Separable, Separable_States, Separable_ExogStates}


struct InitializationOptions{nChoiceVarsZero,C0<:Int}
	"""Objective function at t=0 to choose initial endogenous state variables."""
	problem::Function

	"""Function that maps initial shocks and initial policy into state variables."""
	func::Function

	"""Pointer towards choice of initial state variables."""
	tChoiceVectorsZero::NTuple{nChoiceVarsZero, C0}

	function InitializationOptions{nChoiceVarsZero,C0}(problem::Function,
		func::Function, tChoiceVectorsZero::NTuple{nChoiceVarsZero,C0}) where
			{nChoiceVarsZero,C0<:Int}

		if problem == nothing
			return nothing
		else
			tChoiceVectorsZero[1] == 1 || error(
			"Bad tChoiceVectorsZero: the first state variable must be the (first) choice variable in the intialization problem.")

			new(problem, func, tChoiceVectorsZero)
		end
	end
end
InitializationOptions(problem::Function, func::Function,
	tChoiceVectorsZero::NTuple{nChoiceVarsZero,C0}) where {nChoiceVarsZero,C0} =
	InitializationOptions{nChoiceVarsZero,C0}(problem, func, tChoiceVectorsZero)


struct DiscreteDynamicProblemOptions{RFP<:FuncOrNothing}

	"""The partial reward function."""
	reward_partial::RFP

	"""Options for exact initialization."""
	initialize::Union{InitializationOptions, Nothing}

end
const DDPOptions = DiscreteDynamicProblemOptions


"""
$(TYPEDEF)

Defines an infinite-horizon discrete choice dynamic optimization problem.

# Fields

$(FIELDS)
"""
struct DiscreteDynamicProblem{nStateVars,nChoiceVars,typeC}

	tStateVectors::NTuple{nStateVars, Vector{Float64}} #where N # can use NTuple{N, Vector{Float64}} where N
    tChoiceVectors::NTuple{nChoiceVars, typeC}

    """The reward function defining current period rewards as a function of states and choices."""
    rewardfunc::Function

	"""The transition function of the state variables, depending on the
	integration dimension as a function of states, choices and shocks."""
    transfunc::Function

	"""The integration dimension of the transition function."""
    intdim::Type{ID} where ID<:IntDim # not really necessary, can just have it as a parameter in type...

    """The distribution of the shock(s)."""
    shockdist::Distribution

	"""The discount factor for future rewards."""
    β::Float64 # important to specify here for type stability

	"""Optional information for problem definition."""
	options::Union{DDPOptions, Nothing}


	function DiscreteDynamicProblem{nStateVars,nChoiceVars,typeC}(
		tStateVectors::NTuple{nStateVars, Vector{Float64}},
		tChoiceVectors::NTuple{nChoiceVars, typeC},
		rewardfunc::Function,
		transfunc::Function,
		intdim::Type{<:IntDim},
		shockdist::Distribution,
		β::Float64,
		options::Union{DDPOptions, Nothing}) where
			{nStateVars,nChoiceVars,typeC}


		if intdim <: Separable_Union
			typeC <: Integer || error(
			"Provide a tuple of integers pointing towards the state variables instead of redefining the choice variables.")

			# solver only supports correct order of endogenous state variables
			tChoiceVectors[1] == 1 || error(
			"Bad tChoiceVectors: the first state variable must be the (first) choice variable.")

			if nChoiceVars==2
				tChoiceVectors[2] == 2 || error(
				"Bad tChoiceVectors: the second state variable must be the second choice variable.")
			end
		end

		new{nStateVars,nChoiceVars,typeC}(
			tStateVectors, tChoiceVectors,
			rewardfunc, transfunc, intdim,
			shockdist, β, options)
	end

end
const DDP = DiscreteDynamicProblem
DDP(
	tStateVectors::NTuple{nStateVars, Vector{Float64}},
	tChoiceVectors::NTuple{nChoiceVars, typeC}, args...) where
		{nStateVars,nChoiceVars,typeC} =
	DDP{nStateVars,nChoiceVars,typeC}(
		tStateVectors,
		tChoiceVectors, args...)


"""Options as keyword arguments for convenience."""
function DDP(
			tStateVectors::NTuple{dimStates, Vector{Float64}},
            tChoiceVectors::NTuple{dimChoices, typeC},
            rewardfunc::Function,
            transfunc::Function,
            shockdist::Distribution,
			β::Real;
            intdim::Symbol = :All,
            reward_partial::Union{Function,Nothing} = nothing,
            initializationproblem::Union{Function,Nothing} = nothing,
            initializefunc::Union{Function,Nothing} = nothing,
            tChoiceVectorsZero::Union{NTuple{C0,Int64},Nothing} = nothing,
            ) where {I <: IntDim, dimStates, dimChoices, typeC<:Union{Vector{Float64}, Int64},
                C0}

	intdim in (:All, :Separable, :Separable_ExogStates, :Separable_States) ||
		error("Provided wrong integration dimension $intdim.")

	init_options = InitializationOptions(initializationproblem, initializefunc, tChoiceVectorsZero)

	ddp_options = DDPOptions(reward_partial, init_options)

    DDP(tStateVectors,
        tChoiceVectors,
        rewardfunc,
        transfunc,
        eval(intdim),
        shockdist,
		β,
        ddp_options
        )
end

"""Check whether the integration dimension is separable."""
isseparable(p::DDP) = p.intdim <: Separable_Union

"""Retreive the state vectors that are not choice vectors from the state vector tuple."""
getnonchoicevars(p::DDP) = getnonchoicevars(p.tStateVectors, p.tChoiceVectors)
function getnonchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC}
	allvars = tuple(1:NS...)
	nonchoicevars = tuple(setdiff(Set(allvars), Set(tchoicevars))...)
	return getindex(tStateVectors, collect(nonchoicevars))
end
getnonchoicevars(p::DDP{dimS,dimC,C}) where {dimS,dimC,C<:AbstractVector} = tuple() # i.e. intdim = :All
getnonchoicevars(tStateVectors, tChoiceVectors::Nothing) = nothing
getnonchoicevarszero(p::DDP) = getnonchoicevars(p.tStateVectors, p.options.initialize.tChoiceVectorsZero)


"""Retreive the choice vectors from the state vector tuple."""
getchoicevars(p::DDP) = getchoicevars(p.tStateVectors, p.tChoiceVectors)
getchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC} =
	getindex(tStateVectors, collect(tchoicevars))
"""Retreive nothing if the choice vectors are provided."""
getchoicevars(tStateVectors::NTuple{N1,T1}, tChoiceVectors::NTuple{N2,Vector{T2}}) where
	{N1, T1, N2, T2} = tChoiceVectors
getchoicevars(tStateVectors, tChoiceVectors::Nothing) = nothing
"""Retreive the choice vectors to find policy at t=0."""
getchoicevarszero(p::DDP) = getchoicevars(p.tStateVectors, p.options.initialize.tChoiceVectorsZero)
