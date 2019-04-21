
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

# expand p structure
_simulate(p::DDP{NS,1}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero) where NS =
		_simulate1(sol, shocks, transfunc, p.intdim,
				p.tStateVectors, getchoicevars(p), getnonchoicevars(p),
                getchoicevarszero(p),
				getnonchoicevarszero(p),
				p.options.initialize.func, initialize_exact)


# perhaps better to return a tuple of arrays rather than an array
function _simulate1(sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function,
				intdim::Type{ID},
				tStateVectors::NTuple{dimStates, AbstractVector{T}},
				tChoiceVectors,
				tExogStateVectors::NTuple{dimExogStates, AbstractVector{T}},
				# need information on dimension incase intdim=ExogStates and have only one exog varible
				tChoiceVectorsZero::NTuple{dimChoices0, AbstractVector{T}},
				tExogStateVectorsZero::NTuple{dimExogStates0, AbstractVector{T}},
				initializefunc::Function, initialize_exact::Bool) where
				{T, ID<:IntDim, dimStates, dimExogStates, dimChoices0, dimExogStates0}

	!(initialize_exact && !(typeof(sol) <: DDPSolutionZero)) || error(
		"Solution does not contain intial policy -> rerun solver with `initialize_exact = true`")
	nChoiceVars = 1

	# choose firms as last dimension because loop over firms (easier to code)
	dimShocks, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(dimStates, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)
	# aChoice = fill!(zeros(nChoiceVars, nPeriods+1, nFirms), NaN)

	#############
	# construct interpolator
	# if knots are not evenly spaced: need gridded interpolation
	itp_policy = interpolate(tStateVectors, sol.tmeshPolFun[1], Gridded(Linear()))
	itp_value = interpolate(tStateVectors, sol.meshValFun, Gridded(Linear()))
	itp_policy = extrapolate(itp_policy, Flat())
	itp_value = extrapolate(itp_value, Flat())

	if initialize_exact
        itp_policy0 = interpolate(tExogStateVectorsZero, sol.tmeshPolFunZero[1], Gridded(Linear()))
        itp_policy0 = extrapolate(itp_policy0, Flat())

		function initialize_choices(vExogStateVars0, itp_policy0)
		    C0 = itp_policy0(vExogStateVars0...)
			# C0 = inbounds(C0, tChoiceVectorsZero[1][1], tChoiceVectorsZero[1][end])
        end

        lowerbounds0::NTuple{dimExogStates0,T} = minimum.(tExogStateVectorsZero)
        upperbounds0::NTuple{dimExogStates0,T} = maximum.(tExogStateVectorsZero)
    end

    # get transition of firm i at period t (for period t + 1)
    function fill_it!(vSim_it1, vSim_it, intdim::Type{All}, vChoice, mShocks_it)
        vSim_it1 .= transfunc(vSim_it, vChoice, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable}, vChoice, mShocks_it)
        vSim_it1[1] = vChoice
        vSim_it1[2:end] .= transfunc(vSim_it, vChoice, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_States}, vChoice, mShocks_it)
        vSim_it1[1] = vChoice
        vSim_it1[2:end] .= transfunc(vSim_it, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_ExogStates}, vChoice, mShocks_it)
        vSim_it1[1] = vChoice
        if dimExogStates == 1
            vSim_it1[2] = transfunc(vSim_it[2], mShocks_it) # splat because they're tuples
        else
            vSim_it1[2:end] .= transfunc(vSim_it[2:end], mShocks_it)
        end
    end

    function initialize!(vSim_i1, shocksInit_i)
        if initialize_exact
			if dimShocks > 1
                vSim_i1[2:end] .= inbounds.(initializefunc(shocksInit_i[1]), lowerbounds0, upperbounds0)
			else
                vSim_i1[2:end] .= inbounds.(initializefunc(shocksInit_i), lowerbounds0, upperbounds0)
			end
            vSim_i1[1] = initialize_choices(vSim_i1[2:end], itp_policy0)

        else
            vSim_i1[:] .= initialize_simple(tStateVectors)

        end
    end

	for i = 1:nFirms

		vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[:,:,i]
        mShocks_i = @view shocks.aSim[:,:,i]

        initialize!(@view(vSim_i[:,1]), @view(mShocks_i[:,1]))

		for t = 1 : nPeriods
            # don't need value for now
            # vVal_i[t] = itp_value(vSim_i[:,t]...)

			vChoice = itp_policy(vSim_i[:,t]...)
            # Flat() option in itp_policy ensures it does not leave bounds

			if dimShocks > 1
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t])
			else
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t][1])
			end

			# if do_exog_exit
			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end

	end # firm loop

	return aSim
	# return cat(mVal, aSim, dims=1) #sdata(aSim) # convert shared array to standard array
end # simulatemodel


_simulate(p::DDP{NS,2}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero) where NS =
		_simulate2(sol, shocks, transfunc, p.intdim,
				p.tStateVectors,
                getchoicevars(p), getnonchoicevars(p),
                getchoicevarszero(p), getnonchoicevarszero(p),
				p.options.initialize.func, initialize_exact)
# perhaps better to return a tuple of arrays rather than an array
function _simulate2(sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function,
				intdim::Type{ID},
				tStateVectors::NTuple{dimStates, AbstractVector{T}},
				tChoiceVectors::NTuple{dimChoices, AbstractVector{T}},
				tExogStateVectors::NTuple{dimExogStates, AbstractVector{T}},
                tChoiceVectorsZero::NTuple{dimChoices0, AbstractVector{T}},
				tExogStateVectorsZero::NTuple{dimExogStates0, AbstractVector{T}},
				initializefunc::Function, initialize_exact::Bool) where
					{T, ID<:IntDim, dimStates, dimChoices, dimExogStates, dimChoices0, dimExogStates0}

	!(initialize_exact && !(typeof(sol) <: DDPSolutionZero)) || error(
		"Solution does not contain intial policy -> rerun solver with `initialize_exact = true`")

	# dimChoices::Int64 = length(sol.tmeshPolFun)

	# choose firms as last dimension because loop over firms (easier to code)
	dimShocks, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(dimStates, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)

	#############
	# construct interpolator
	# if knots are not evenly spaced: need gridded interpolation
	itp_policy = [interpolate(tStateVectors, polfun, Gridded(Linear())) for polfun in sol.tmeshPolFun]
	itp_value = interpolate(tStateVectors, sol.meshValFun, Gridded(Linear()))

	# if knots are evenly spaced: scaled Bsplines (not working yet)
	# itp = interpolate(meshPolFun, BSpline(Linear()), OnGrid())
	# itp = scale(itp, vK, vAlog)

	# specify that if value requested outside of grid, return grid value
	# I am using my inbounds() function, but still there is an issue with values on the grid somehow..
	# hope to remove these lines in the future
	itp_policy = [extrapolate(itp, Flat()) for itp in itp_policy]
	itp_value = extrapolate(itp_value, Flat())

	if initialize_exact
        itp_policy0 = [interpolate(tExogStateVectorsZero, polfun, Gridded(Linear())) for polfun in sol.tmeshPolFunZero]
        itp_policy0 = [extrapolate(itp, Flat()) for itp in itp_policy0]

        initialize_choices(vExogStateVars0) = [itp(vExogStateVars0...) for itp in itp_policy0]

        lowerbounds0::NTuple{dimExogStates0,T} = minimum.(tExogStateVectorsZero)
        upperbounds0::NTuple{dimExogStates0,T} = maximum.(tExogStateVectorsZero)
    end

	# # for exogenous exit
 	# do_exog_exit = isdefined(params,:π)
	# if do_exog_exit
	# 	do_exog_exit =  do_exog_exit * (params.π>0)
	# end

    # get transition of firm i at period t (for period t + 1)
    function fill_it!(vSim_it1, vSim_it, intdim::Type{All}, vChoice, mShocks_it)
        vSim_it1 .= transfunc(vSim_it, vChoice, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable}, vChoice, mShocks_it)
        vSim_it1[1:dimChoices] .= vChoice
        vSim_it1[1+dimChoices:end] .= transfunc(vSim_it, vChoice, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it,intdim::Type{Separable_States}, vChoice, mShocks_it)
        vSim_it1[1:dimChoices] .= vChoice
        vSim_it1[1+dimChoices:end] .= transfunc(vSim_it, mShocks_it)
    end
    function fill_it!(vSim_it1, vSim_it, intdim::Type{Separable_ExogStates}, vChoice, mShocks_it)
        vSim_it1[1:dimChoices] .= vChoice
        if dimExogStates == 1
            vSim_it1[1+dimChoices:end] .= transfunc(vSim_it[1+dimChoices], mShocks_it)
        else
            vSim_it1[1+dimChoices:end] .= transfunc(vSim_it[1+dimChoices:end], mShocks_it)
        end
    end

    function initialize!(vSim_i1, shocksInit_i)
        if initialize_exact
			if dimShocks > 1
                vSim_i1[1+dimChoices0:end] .= inbounds.(initializefunc(shocksInit_i[1]), lowerbounds0, upperbounds0)
			else
                vSim_i1[1+dimChoices0:end] .= inbounds.(initializefunc(shocksInit_i), lowerbounds0, upperbounds0)
			end
            vSim_i1[1:dimChoices0] = initialize_choices(vSim_i1[1+dimChoices0:end])

        else
            vSim_i1[:] .= initialize_simple(tStateVectors) # stable

        end
    end


	for i = 1:nFirms

        vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[:,:,i]
        mShocks_i = @view shocks.aSim[:,:,i]

        initialize!(@view(vSim_i[:,1]), @view(mShocks_i[:,1]))

		for t = 1 : nPeriods
            # don't need value for now
            # vVal_i[t] = itp_value(vSim_i[:,t]...)

			vChoice = [itp(vSim_i[:,t]...) for itp in itp_policy]

			if dimShocks > 1
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t])
			else
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t][1])
			end

			# if do_exog_exit
			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end

	end # firm loop

	return aSim#, itp_policy
	# return cat(mVal, aSim, dims=1) #sdata(aSim) # convert shared array to standard array
end # simulatemodel
