# convenience wrappers
simulate(sol::AbstractDDPSolution;
			nPeriods::Int64 = 60,
			nFirms::Int64 = 100, kwargs...) =
	_simulate(sol.prob, sol, drawshocks(sol.prob, nPeriods=nPeriods, nFirms=nFirms),
				sol.prob.transfunc; kwargs...)
simulate(sol::AbstractDDPSolution, shocks::DDPShocks; kwargs...) =
	_simulate(sol.prob, sol, shocks, sol.prob.transfunc; kwargs...)

simulate(p::DDP; kwargs...) =
	_simulate(p, drawshocks(p); kwargs...)
simulate(p::DDP, shocks::DDPShocks; kwargs...) =
	_simulate(p, solve(p), shocks, p.transfunc; kwargs...)


function initialize_simple(tStateVectors::NTuple{N, Vector{Float64}}) where N
	# write a standard initialization function if want burn-in method
	# get the middle index of all state variables to initialize
	getindex.(tStateVectors, Int.(floor.(length.(tStateVectors)./2)))
end

# function fill_it!(p::DDP{NS,NC,C,All}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,NC,C}
# 	vSim_it1 .= p.transfunc(vSim_it, vChoice, mShocks_it)
# end

# function fill_it!(p::DDP{NS,1,C,Separable}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,C}
# 	vSim_it1[1] = vChoice
# 	vSim_it1[2:end] .= p.transfunc(vSim_it, vChoice, mShocks_it)
# end
# function fill_it!(p::DDP{NS,NC,C,Separable}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,NC,C}
# 	vSim_it1[1:NC] .= vChoice
# 	vSim_it1[1+NC:end] .= p.transfunc(vSim_it, vChoice, mShocks_it)
# end

# function fill_it!(p::DDP{NS,1,C,Separable_States}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,C}
# 	vSim_it1[1] = vChoice
# 	vSim_it1[2:end] .= p.transfunc(vSim_it, mShocks_it)
# end
# function fill_it!(p::DDP{NS,NC,C,Separable_States}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,NC,C}
# 	vSim_it1[1:NC] .= vChoice
# 	vSim_it1[1+NC:end] .= p.transfunc(vSim_it, mShocks_it)
# end

# function fill_it!(p::DDP{NS,1,C,Separable_ExogStates}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,C}
# 	vSim_it1[1] = vChoice
# 	# if dimExogStates == 1
# 		vSim_it1[2] = p.transfunc(vSim_it[2], mShocks_it) # splat because they're tuples
# 	# else
# 		# vSim_it1[2:end] .= p.transfunc(vSim_it[2:end], mShocks_it)
# 	# end
# end
# function fill_it!(p::DDP{NS,NC,C,Separable_ExogStates}, vSim_it1, vSim_it, vChoice, mShocks_it) where {NS,NC,C}
# 	vSim_it1[1:NC] .= vChoice
# 	# if dimExogStates == 1
# 		vSim_it1[1+NC:end] .= p.transfunc(vSim_it[1+NC], mShocks_it)
# 	# else
# 	# 	vSim_it1[1+NC:end] .= transfunc(vSim_it[1+NC:end], mShocks_it)
# 	# end
# end


# expand p structure
_simulate(p::DDP{NS,1}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero,
		get_value::Bool = false) where NS =
		_simulate(p, sol, shocks, transfunc, p.intdim,
				p.tStateVectors, getchoicevars(p), getnonchoicevars(p),
                getchoicevarszero(p),
				getnonchoicevarszero(p),
				p.options.initialize.func,
				initialize_exact,
				get_value)


# perhaps better to return a tuple of arrays rather than an array
function _simulate(p::DDP{dimStates,1}, sol::AbstractDDPSolution, shocks::DDPShocks{dimShocks}, transfunc,
				intdim::Type{ID},
				tStateVectors::NTuple{dimStates, AbstractVector{T}},
				tChoiceVectors,
				tExogStateVectors::NTuple{dimExogStates, AbstractVector{T}},
				# need information on dimension incase intdim=ExogStates and have only one exog variable
				tChoiceVectorsZero::NTuple{dimChoices0, AbstractVector{T}},
				tExogStateVectorsZero::NTuple{dimExogStates0, AbstractVector{T}},
				initializefunc, initialize_exact::Bool,
				get_value::Bool) where
				{T, ID<:IntDim, dimStates, dimExogStates, dimChoices0,
					dimExogStates0, dimShocks}

	!(initialize_exact && !(typeof(sol) <: DDPSolutionZero)) || error(
		"Solution does not contain intial policy -> rerun solver with `initialize_exact = true`")
	nChoiceVars = 1

	# choose firms as last dimension because loop over firms (easier to code)
	unused, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(dimStates, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)
	aChoice = fill!(zeros(nChoiceVars, nPeriods+1, nFirms), NaN)

	value = sol.value
	policy = sol.policy[1]

	if initialize_exact
		policy0 = sol.policy0[1]

		initialize_choices(vExogStateVars0, policy0) = policy0(vExogStateVars0...)

        lowerbounds0::NTuple{dimExogStates0,T} = minimum.(tExogStateVectorsZero)
        upperbounds0::NTuple{dimExogStates0,T} = maximum.(tExogStateVectorsZero)
    end

    # get transition of firm i at period t (for period t + 1)
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{All}, vChoice, mShocks_it)
    #     vSim_it1 .= transfunc(vSim_it, vChoice, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable}, vChoice, mShocks_it)
    #     vSim_it1[1] = vChoice
    #     vSim_it1[2:end] .= transfunc(vSim_it, vChoice, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_States}, vChoice, mShocks_it)
    #     vSim_it1[1] = vChoice
    #     vSim_it1[2:end] .= transfunc(vSim_it, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_ExogStates}, vChoice, mShocks_it)
    #     vSim_it1[1] = vChoice
    #     if dimExogStates == 1
    #         vSim_it1[2] = transfunc(vSim_it[2], mShocks_it) # splat because they're tuples
    #     else
    #         vSim_it1[2:end] .= transfunc(vSim_it[2:end], mShocks_it)
    #     end
    # end

    function initialize!(vSim_i1, shocksInit_i)
        if initialize_exact
			if dimShocks == 1
                # vSim_i1[2:end] .= inbounds.(initializefunc(shocksInit_i[1]), lowerbounds0, upperbounds0)
                vSim_i1[2:end] .= initializefunc(shocksInit_i[1])
			else
                # vSim_i1[2:end] .= inbounds.(initializefunc(shocksInit_i), lowerbounds0, upperbounds0)
				# @show shocksInit_i
				vSim_i1[2:end] .= initializefunc(shocksInit_i)
			end
            vSim_i1[1] = initialize_choices(vSim_i1[2:end], policy0)
        else
            vSim_i1[:] .= initialize_simple(tStateVectors)
        end
    end

	for i = 1:nFirms

		vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[:,:,i]
		vChoice_i = @view aChoice[:,:,i]
        mShocks_i = @view shocks.aSim[:,:,i]

        initialize!(@view(vSim_i[:,1]), @view(mShocks_i[:,1]))

		for t = 1 : nPeriods

            if get_value
            	vVal_i[t] = value(vSim_i[:,t]...)
			end

			vChoice_i[t] = policy(vSim_i[:,t]...)

			if dimShocks > 1
				get_transition!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice_i[t], mShocks_i[:,t])
			else
				# fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice_i[t], mShocks_i[:,t][1])
				# fill_it!(p, @view(vSim_i[:,t+1]), vSim_i[:,t], vChoice_i[t], mShocks_i[:,t][1])
				get_transition!(p.transfunc, @view(vSim_i[:,t+1]), vSim_i[:,t], vChoice_i[t], mShocks_i[:,t][1])
			end

			# if do_exog_exit
			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end

	end # firm loop

	if get_value
		return DDPSimulation(p, sol, mVal, aSim, aChoice)
	else
		return DDPSimulation(p, sol, nothing, aSim, aChoice)
	end
end # simulatemodel


_simulate(p::DDP{NS,2}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero,
		get_value::Bool = false) where NS =
		_simulate(p, sol, shocks, transfunc, p.intdim,
				p.tStateVectors,
                getchoicevars(p), getnonchoicevars(p),
                getchoicevarszero(p), getnonchoicevarszero(p),
				p.options.initialize.func, initialize_exact,
				get_value)
# perhaps better to return a tuple of arrays rather than an array
function _simulate(p::DDP{dimStates,2}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc,
				intdim::Type{ID},
				tStateVectors::NTuple{dimStates, AbstractVector{T}},
				tChoiceVectors::NTuple{dimChoices, AbstractVector{T}},
				tExogStateVectors::NTuple{dimExogStates, AbstractVector{T}},
                tChoiceVectorsZero::NTuple{dimChoices0, AbstractVector{T}},
				tExogStateVectorsZero::NTuple{dimExogStates0, AbstractVector{T}},
				initializefunc, initialize_exact::Bool,
				get_value::Bool) where
					{T, ID<:IntDim, dimStates, dimChoices, dimExogStates, dimChoices0, dimExogStates0}

	!(initialize_exact && !(typeof(sol) <: DDPSolutionZero)) || error(
		"Solution does not contain intial policy -> rerun solver with `initialize_exact = true`")

	# choose firms as last dimension because loop over firms (easier to code)
	dimShocks::Int, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(dimStates, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)
	aChoice = fill!(zeros(dimChoices, nPeriods+1, nFirms), NaN)

	value = sol.value
	policy = sol.policy

	if initialize_exact
        policy0 = sol.policy0

        initialize_choices(vExogStateVars0) = [itp(vExogStateVars0...) for itp in policy0]

        lowerbounds0::NTuple{dimExogStates0,T} = minimum.(tExogStateVectorsZero)
        upperbounds0::NTuple{dimExogStates0,T} = maximum.(tExogStateVectorsZero)
    end

	# # for exogenous exit
 	# do_exog_exit = isdefined(params,:π)
	# if do_exog_exit
	# 	do_exog_exit =  do_exog_exit * (params.π>0)
	# end

    # get transition of firm i at period t (for period t + 1)
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{All}, vChoice, mShocks_it)
    #     vSim_it1 .= transfunc(vSim_it, vChoice, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable}, vChoice, mShocks_it)
    #     vSim_it1[1:dimChoices] .= vChoice
    #     vSim_it1[1+dimChoices:end] .= transfunc(vSim_it, vChoice, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it,intdim::Type{Separable_States}, vChoice, mShocks_it)
    #     vSim_it1[1:dimChoices] .= vChoice
    #     vSim_it1[1+dimChoices:end] .= transfunc(vSim_it, mShocks_it)
    # end
    # function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_ExogStates}, vChoice, mShocks_it)
    #     vSim_it1[1:dimChoices] .= vChoice
    #     if dimExogStates == 1
    #         vSim_it1[1+dimChoices:end] .= transfunc(vSim_it[1+dimChoices], mShocks_it)
    #     else
    #         vSim_it1[1+dimChoices:end] .= transfunc(vSim_it[1+dimChoices:end], mShocks_it)
    #     end
    # end

    function initialize!(vSim_i1, shocksInit_i)
        if initialize_exact
			if dimShocks == 1
                # vSim_i1[1+dimChoices0:end] .= inbounds.(initializefunc(shocksInit_i[1]), lowerbounds0, upperbounds0)
                vSim_i1[1+dimChoices0:end] .= initializefunc(shocksInit_i[1])
			else
                # vSim_i1[1+dimChoices0:end] .= inbounds.(initializefunc(shocksInit_i), lowerbounds0, upperbounds0)
                vSim_i1[1+dimChoices0:end] .= initializefunc(shocksInit_i)
			end
            vSim_i1[1:dimChoices0] = initialize_choices(vSim_i1[1+dimChoices0:end])

        else
            vSim_i1[:] .= initialize_simple(tStateVectors) # stable

        end
    end


	for i = 1:nFirms

        vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[:,:,i]
		vChoice_i = @view aChoice[:,:,i]
        mShocks_i = @view shocks.aSim[:,:,i]

        initialize!(@view(vSim_i[:,1]), @view(mShocks_i[:,1]))

		for t = 1 : nPeriods
			if get_value
            	vVal_i[t] = value(vSim_i[:,t]...)
			end

			for j = 1 : dimChoices
				vChoice_i[j,t] = policy[j](vSim_i[:,t]...)
			end

			if dimShocks > 1
				get_transition!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice_i[:,t], mShocks_i[:,t])
			else
				# fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice_i[:,t], mShocks_i[:,t][1])
				# fill_it!(p, @view(vSim_i[:,t+1]), vSim_i[:,t], vChoice_i[:,t], mShocks_i[:,t][1])
				get_transition!(p.transfunc, @view(vSim_i[:,t+1]), vSim_i[:,t], vChoice_i[:,t], mShocks_i[:,t][1])
			end

			# if do_exog_exit
			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end

	end # firm loop

	if get_value
		return DDPSimulation(p, sol, mVal, aSim, aChoice)
	else
		return DDPSimulation(p, sol, nothing, aSim, aChoice)
	end
end # simulatemodel
