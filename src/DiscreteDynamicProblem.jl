
# model definition and parameters
# abstract type DiscreteDynamicModel end
# const DDP = DiscreteDynamicModel

struct InitializationOptions{NC0,NSh,C0<:Int,D<:Distribution,P,F}
	"""Objective function at t=0 to choose initial endogenous state variables."""
	problem::P

	"""Function that maps initial shocks and initial policy into state variables."""
	func::F

	"""Distribution of initial shocks."""
	shockdist::D

	"""Pointer towards choice of initial state variables."""
	tChoiceVectorsZero::NTuple{NC0, C0}

	function InitializationOptions{NC0,NSh,C0,D,P,F}(problem,
		func, shockdist,
		tChoiceVectorsZero::NTuple{NC0,C0},
		tStateVectors::NTuple{NS,Vector{Float64}}) where
			{NC0,NSh,C0<:Int,NS,D,P,F}

		tChoiceVectorsZero[1] == 1 || error(
		"Bad tChoiceVectorsZero: the first state variable must be the (first) choice variable in the intialization problem.")

		tExogStateVectorsZero = getnonchoicevars(tStateVectors, tChoiceVectorsZero)
		func_inbounds = wrapinbounds(func, tExogStateVectorsZero)

		new{NC0,NSh,C0,typeof(shockdist),typeof(problem),typeof(func_inbounds)}(
			problem, func_inbounds, shockdist, tChoiceVectorsZero)

	end
end
InitializationOptions(problem, func, shockdist::Distribution,
	tChoiceVectorsZero::NTuple{NC0,C0},
	tStateVectors::NTuple{NS,Vector{Float64}}) where {NC0,C0,NS} =
	InitializationOptions{NC0,length(shockdist),C0,typeof(shockdist),
		typeof(problem),typeof(func)}(
		problem, func, shockdist, tChoiceVectorsZero,tStateVectors)


struct DiscreteDynamicProblemOptions{RFP, I<:Union{InitializationOptions, Nothing}}

	"""The partial reward function."""
	rewardfunc_partial::RFP

	"""Options for exact initialization."""
	initialize::I

end
const DDPOptions = DiscreteDynamicProblemOptions


"""
$(TYPEDEF)

Defines an infinite-horizon discrete choice dynamic optimization problem.

# Fields

$(FIELDS)
"""
struct DiscreteDynamicProblem{nStateVars,nChoiceVars,typeC,ID<:IntDim,RF,TF,D<:Distribution}

	tStateVectors::NTuple{nStateVars, Vector{Float64}}
    tChoiceVectors::NTuple{nChoiceVars, typeC}

    """The reward function defining current period rewards as a function of states and choices."""
    rewardfunc::RF

	"""The transition function of the state variables, depending on the
	integration dimension as a function of states, choices and shocks."""
    transfunc::TF

	"""The integration dimension of the transition function."""
    intdim::Type{ID}

    """The distribution of the shock(s)."""
    shockdist::D

	"""The discount factor for future rewards."""
    β::Float64

	"""Optional information for problem definition."""
	options::DDPOptions


	function DiscreteDynamicProblem{nStateVars,nChoiceVars,typeC}(
		tStateVectors::NTuple{nStateVars, Vector{Float64}},
		tChoiceVectors::NTuple{nChoiceVars, typeC},
		rewardfunc,
		transfunc,
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


		# extend transitionfunc to stay inbounds
		# note: not sure this will yield perfect performance
		# perhaps better to do the wrapping where the transition function is used
		if intdim == All
			transfunc_inbounds = wrapinbounds(transfunc, tStateVectors)
			nExogStates = 0
	    else
			tExogStateVectors = getnonchoicevars(tStateVectors, tChoiceVectors)
			transfunc_inbounds = wrapinbounds(transfunc, tExogStateVectors)
			nExogStates = length(tExogStateVectors)
	    end

		nShocks = length(shockdist)
		transition = Transition{intdim,nChoiceVars,nShocks,
			nExogStates,typeof(transfunc_inbounds)}(transfunc_inbounds)

		new{nStateVars,nChoiceVars,typeC,intdim,typeof(rewardfunc),typeof(transition),typeof(shockdist)}(
			tStateVectors, tChoiceVectors,
			rewardfunc, transition, intdim,
			shockdist, β, options)
		# new{nStateVars,nChoiceVars,typeC,intdim,typeof(rewardfunc),typeof(transfunc_inbounds),typeof(shockdist)}(
		# 	tStateVectors, tChoiceVectors,
		# 	rewardfunc, transfunc_inbounds, intdim,
		# 	shockdist, β, options)
	end

end



const DDP = DiscreteDynamicProblem
DDP(tStateVectors::NTuple{nStateVars, Vector{Float64}},
	tChoiceVectors::NTuple{nChoiceVars, typeC}, args...) where
		{nStateVars,nChoiceVars,typeC} =
	DDP{nStateVars,nChoiceVars,typeC}(
		tStateVectors,
		tChoiceVectors, args...)


"""Options as keyword arguments for convenience."""
function DDP(
			tStateVectors::NTuple{dimStates, Vector{Float64}},
            tChoiceVectors::NTuple{dimChoices, typeC},
            rewardfunc,
            transfunc,
            shockdist::Distribution,
			β::Real;
            intdim::Symbol = :All,
            rewardfunc_partial = nothing,
            initializationproblem = nothing,
            initializefunc = nothing,
			shockdist_initial::Union{Distribution,Nothing} = shockdist,
            tChoiceVectorsZero::Union{NTuple{C0,Int64},Nothing} = nothing,
            ) where {I <: IntDim, dimStates, dimChoices, typeC<:Union{Vector{Float64}, Int64},
                C0}

	# make sure functions are callable
	!isempty(methods(rewardfunc)) || error("supplied rewardfunc is not callable")
	!isempty(methods(transfunc)) || error("supplied transfunc is not callable")

	intdim in (:All, :Separable, :Separable_ExogStates, :Separable_States) ||
		error("Provided wrong integration dimension $intdim.")

	if initializationproblem == nothing
		init_options = nothing
	else
		# tExogStateVectorsZero = getnonchoicevars(tStateVectors, tChoiceVectorsZero)
		# initializefunc_inbounds = wrapinbounds(initializefunc, tExogStateVectorsZero)
		# init_options = InitializationOptions(initializationproblem, initializefunc_inbounds,
		# 	tChoiceVectorsZero)
		init_options = InitializationOptions(initializationproblem, initializefunc,
			shockdist_initial, tChoiceVectorsZero, tStateVectors)
	end

	ddp_options = DDPOptions(rewardfunc_partial, init_options)

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
"""Retreive exogenous states at t=0."""
function getnonchoicevarszero(p::DDP)
	if p.options.initialize == nothing
		return nothing
	else
		return 	getnonchoicevars(p.tStateVectors, p.options.initialize.tChoiceVectorsZero)
	end
end

"""Retreive the choice vectors from the state vector tuple."""
getchoicevars(p::DDP) = getchoicevars(p.tStateVectors, p.tChoiceVectors)
getchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC} =
	getindex(tStateVectors, collect(tchoicevars))
"""Retreive nothing if the choice vectors are provided."""
getchoicevars(tStateVectors::NTuple{N1,T1}, tChoiceVectors::NTuple{N2,Vector{T2}}) where
	{N1, T1, N2, T2} = tChoiceVectors
getchoicevars(tStateVectors, tChoiceVectors::Nothing) = nothing
"""Retreive endogenous states at t=0."""
function getchoicevarszero(p::DDP)
	if p.options.initialize == nothing
		return nothing
	else
		return getchoicevars(p.tStateVectors, p.options.initialize.tChoiceVectorsZero)
	end
end
