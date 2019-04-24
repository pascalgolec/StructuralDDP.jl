
"""Initialize the simulation from the middle of the state space."""
initialize_simple!(vStates0, sol::AbstractDDPSolution) =
	initialize_simple!(vStates0, sol.prob.tStateVectors)
function initialize_simple!(vStates0, tStateVectors::NTuple{N, Vector{Float64}}) where N
	vStates0 .= getindex.(tStateVectors, Int.(floor.(length.(tStateVectors)./2)))
end

initialize_choices(sol::DDPSolutionZero{NS,NC,NE0,NC0}, vExogStates0) where {NS,NC,NE0,NC0} =
	[itp(vExogStates0...) for itp in sol.policy0]
initialize_choices(sol::DDPSolutionZero{NS,NC,NE0,1}, vExogStates0) where {NS,NC,NE0} =
	sol.policy0[1](vExogStates0...)

"""Initialization with shocks and choices"""
function initialize_exact!(vStates0, sol::DDPSolutionZero{NS,NC,NE0,NC0},
		init_opt::InitializationOptions{NC0}, vShocks) where {NS,NC,NE0,NC0}
	vStates0[1+NC0:end] .= init_opt.func(vShocks)
	vStates0[1:NC0] .= initialize_choices(sol, vStates0[1+NC0:end])
end
"""when initial shock is one-dimensional."""
function initialize_exact!(vStates0, sol::DDPSolutionZero{NS,NC,NE0,NC0},
		init_opt::InitializationOptions{NC0,1}, vShocks) where {NS,NC,NE0,NC0}
	vStates0[1+NC0:end] .= init_opt.func(vShocks[1])
	vStates0[1:NC0] .= initialize_choices(sol, vStates0[1+NC0:end])
end
"""Only one choice variable"""
function initialize_exact!(vStates0, sol::DDPSolutionZero{NS,NC,NE0,1},
		init_opt::InitializationOptions{1}, vShocks) where {NS,NC,NE0}
	vStates0[2:end] .= init_opt.func(vShocks)
	vStates0[1] = initialize_choices(sol, vStates0[2:end])
end
"""Only one choice variable and shock is one-dimensional."""
function initialize_exact!(vStates0, sol::DDPSolutionZero{NS,NC,NE0,1},
		init_opt::InitializationOptions{1,1}, vShocks) where {NS,NC,NE0}
	vStates0[2:end] .= init_opt.func(vShocks[1])
	vStates0[1] = initialize_choices(sol, vStates0[2:end])
end


function initialize!(vStates0, initialize_exact::Bool, sol::AbstractDDPSolution,
		init_opt, vShocks)
	if initialize_exact
		initialize_exact!(vStates0, sol, init_opt, vShocks)
	else
		initialize_simple!(vStates0, sol)
	end
end


"""Get optimal choice of a firm at time t."""
function get_choice_t!(vChoice, policy::NTuple{NC,T}, vStates, t) where {NC,T}
	for j = 1 : NC
		vChoice[j,t] = policy[j](vStates[:,t]...)
	end
end

"""Get transition for a firm at time t."""
get_transition_t!(transfunc::Transition,
	vStates, vChoices, vShocks, t) =
	get_transition!(transfunc, @view(vStates[:,t+1]), vStates[:,t],
		vChoices[:,t], vShocks[:,t])
"""Only one shock."""
get_transition_t!(transfunc::Transition{ID,NC,1},
	vStates, vChoices, vShocks, t) where {ID,NC} =
	get_transition!(transfunc, @view(vStates[:,t+1]), vStates[:,t],
		vChoices[:,t], vShocks[:,t][1])
"""Only one choice var."""
get_transition_t!(transfunc::Transition{ID,1},
	vStates, vChoices, vShocks, t) where ID =
	get_transition!(transfunc, @view(vStates[:,t+1]), vStates[:,t],
		vChoices[t], vShocks[:,t])
"""Only one choice var and only one shock."""
get_transition_t!(transfunc::Transition{ID,1,1},
	vStates, vChoices, vShocks, t) where ID =
	get_transition!(transfunc, @view(vStates[:,t+1]), vStates[:,t],
		vChoices[t], vShocks[:,t][1])

# convenience wrappers
simulate(sol::AbstractDDPSolution;
			nPeriods::Int64 = 60,
			nFirms::Int64 = 100, kwargs...) =
	_simulate(sol.prob, sol, drawshocks(sol.prob, nPeriods=nPeriods, nFirms=nFirms); kwargs...)

simulate(sol::AbstractDDPSolution, shocks::AbstractDDPShocks; kwargs...) =
	_simulate(sol.prob, sol, shocks; kwargs...)

simulate(p::DDP; kwargs...) =
	_simulate(p, drawshocks(p); kwargs...)
simulate(p::DDP, shocks::AbstractDDPShocks; kwargs...) =
	_simulate(p, solve(p), shocks; kwargs...)


_simulate(p::DDP, sol::AbstractDDPSolution, shocks::AbstractDDPShocks;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero,
		get_value::Bool = false) =
		_simulate(p, sol, shocks, p.transfunc, p.options.initialize,
				initialize_exact, get_value)

function _simulate(p::DDP{NS,NC}, sol::AbstractDDPSolution, shocks::AbstractDDPShocks,
				transfunc::Transition{ID,NC,dimShocks},
				init_opt::Union{InitializationOptions,Nothing},
				initialize_exact::Bool,
				get_value::Bool) where {ID, NS, NC, dimShocks}

	!(initialize_exact && !(typeof(sol) <: DDPSolutionZero)) || error(
		"Solution does not contain intial policy -> specify initializationproblem when defining DDP")

	unused, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(NS, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)
	aChoice = fill!(zeros(NC, nPeriods+1, nFirms), NaN)

	value = sol.value
	policy = sol.policy

	# # for exogenous exit
 	# do_exog_exit = isdefined(params,:π)
	# if do_exog_exit
	# 	do_exog_exit =  do_exog_exit * (params.π>0)
	# end

	for i = 1:nFirms

        vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[:,:,i]
		vChoice_i = @view aChoice[:,:,i]
        mShocks_i = @view shocks.aSim[:,:,i]

		if initialize_exact
        	mShocks0_i = @view shocks.aInit[:,1,i]
		else
			mShocks0_i = nothing
		end

		initialize!(@view(vSim_i[:,1]), initialize_exact, sol, init_opt,
				 mShocks0_i)

		for t = 1 : nPeriods

			if get_value
            	vVal_i[t] = value(vSim_i[:,t]...)
			end

			get_choice_t!(vChoice_i, policy, vSim_i, t)

			get_transition_t!(transfunc, vSim_i, vChoice_i, mShocks_i, t)

			# if do_exog_exit
			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end # period t

	end # firm i

	if get_value
		return DDPSimulation(p, sol, mVal, aSim, aChoice)
	else
		return DDPSimulation(p, sol, nothing, aSim, aChoice)
	end
end # simulatemodel
