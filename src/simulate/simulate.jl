
# function simulate(p::DDM)
# 	# seed!(42) # fix the seed
# 	simulate(p, drawshocks(p))
# end

# function simulate(p::DDM, sol::DDPSolution)
# 	# seed!(42) # fix the seed
# 	simulate(p, sol, drawshocks(p))
# end

# convenience wrappers
simulate(p::DDM; kwargs...) =
	_simulate(p, drawshocks(p); kwargs...)
simulate(p::DDM, sol::AbstractDDPSolution;
			nPeriods::Int64 = 60,
			nFirms::Int64 = 100, kwargs...) =
	_simulate(p, sol, drawshocks(p, nPeriods=nPeriods, nFirms=nFirms),
				p.transfunc; kwargs...)
simulate(p::DDM, shocks::DDPShocks; kwargs...) =
	_simulate(p, solve(p), shocks, p.transfunc; kwargs...)
simulate(p::DDM, sol::AbstractDDPSolution, shocks::DDPShocks; kwargs...) =
	_simulate(p, sol, shocks, p.transfunc; kwargs...)

function initialize_simple(tStateVectors::NTuple{N, Vector{Float64}}) where N
	# write a standard initialization function if want burn-in method
	# get the middle index of all state variables to initialize
	# convert(Vector{Float64},
	# 	getindex.(p.tStateVectors, Int.(floor.(length.(p.tStateVectors)./2))))
	getindex.(tStateVectors, Int.(floor.(length.(tStateVectors)./2)))
end

# expand p structure
_simulate(p::DDP{NS,1}, sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function;
        initialize_exact::Bool = typeof(sol) <: DDPSolutionZero) where NS =
		_simulate1(sol, shocks, transfunc, p.intdim,
				p.tStateVectors, getchoicevars(p), getnonchoicevars(p),
				getnonchoicevarszero(p),
				p.initializefunc, initialize_exact)


# perhaps better to return a tuple of arrays rather than an array
function _simulate1(sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function,
				intdim::Type{T},
				tStateVectors::NTuple{dimStates, AbstractVector},
				tChoiceVectors,
				tExogStateVectors::NTuple{dimExogStates, AbstractVector},
				# need information on dimension incase intdim=ExogStates and have only one exog varible
				tOtherStateVectorsZero,
				initializefunc::Function, initialize_exact::Bool) where
				{T<:DDMIntDim, dimStates, dimExogStates}

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
        itp_policy0 = interpolate(tOtherStateVectorsZero, sol.tmeshPolFunZero[1], Gridded(Linear()))
        itp_policy0 = extrapolate(itp_policy0, Flat())
    end

	minChoice::Float64 = minimum(tChoiceVectors[1])
	maxChoice::Float64 = maximum(tChoiceVectors[1])

	for i = 1:nFirms

		vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[1,:,i]

        if initialize_exact
            vSim_i[:,1] .= initializefunc(shocks.aInit[:,1,i], itp_policy0)
        else
            vSim_i[:,1] .= initialize_simple(tStateVectors) # stable
        end

		vVal_i[1] = itp_value(vSim_i[:,1]...)

		mShocks_i = @view shocks.aSim[:,:,i]


		for t = 1 : nPeriods
			vChoice = inbounds(itp_policy(vSim_i[:,t]...), minChoice, maxChoice)

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
					vSim_it1[2:end] .= transfunc(vSim_it[2], mShocks_it)
				else
					vSim_it1[2:end] .= transfunc(vSim_it[2:end], mShocks_it)
				end
			end

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
				p.tStateVectors, getchoicevars(p), getnonchoicevars(p),
				getnonchoicevarszero(p),
				p.initializefunc, initialize_exact)
# perhaps better to return a tuple of arrays rather than an array
function _simulate2(sol::AbstractDDPSolution, shocks::DDPShocks, transfunc::Function,
				intdim::Type{T},
				tStateVectors::NTuple{dimStates, AbstractVector},
				tChoiceVectors::NTuple{dimChoices, AbstractVector},
				tExogStateVectors::NTuple{dimExogStates, AbstractVector},
				tOtherStateVectorsZero,
				initializefunc::Function, initialize_exact::Bool) where
					{T<:DDMIntDim, dimStates, dimChoices, dimExogStates}

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
        itp_policy0 = [interpolate(tOtherStateVectorsZero, polfun, Gridded(Linear())) for polfun in sol.tmeshPolFunZero]
        itp_policy0 = [extrapolate(itp, Flat()) for itp in itp_policy0]
    end

	# # for exogenous exit
 	# do_exog_exit = isdefined(params,:π)
	# if do_exog_exit
	# 	do_exog_exit =  do_exog_exit * (params.π>0)
	# end

	# simulate one firm at a time, so that can do in parallel later
	# @sync @parallel for i = 1:nFirms

	minChoice::Vector{Float64} = [minimum(choicevec) for choicevec in tChoiceVectors]
	maxChoice::Vector{Float64} = [maximum(choicevec) for choicevec in tChoiceVectors]

	for i = 1:nFirms

		vSim_i = @view aSim[:,:,i]
		vVal_i = @view mVal[1,:,i]

		# vSim_i = fill!(zeros(nDim, nPeriods+1), NaN) # need NaN so that will register as exit
		# vVal_i = fill!(zeros(nPeriods+1), NaN)

        if initialize_exact
            vSim_i[:,1] .= initializefunc(shocks.aInit[:,1,i], itp_policy0...)
        else
            vSim_i[:,1] .= initialize_simple(tStateVectors) # stable
        end

		vVal_i[1] = itp_value(vSim_i[:,1]...)

		mShocks_i = @view shocks.aSim[:,:,i]


		for t = 1 : nPeriods

			vChoice = [inbounds(itp_policy[n](vSim_i[:,t]...), minChoice[n], maxChoice[n])
				for n = 1:dimChoices]

			# THIS BELOW LEADS TO TYPE INFERENCE PROBLEMS or is just slow
			# vSim_i[:,t+1] .= transfunc(SA, vSim_i[:,t], choice, mShocks_i[:,t])

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

			if dimShocks > 1
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t])
			else
				fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], intdim, vChoice, mShocks_i[:,t][1])
			end

			# if dimShocks > 1
			# 	fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], transfunc, intdim, vChoice, mShocks_i[:,t], dimChoices)
			# else
			# 	fill_it!(@view(vSim_i[:,t+1]), vSim_i[:,t], transfunc, intdim, vChoice, mShocks_i[:,t][1], dimChoices)
			# end

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

# can be a bit faster if have a separate method for singlechoicevar problems. then need to index less...
# the gains while testing were minimal though
# function _simulate(sol::DDPSolution{1}, shocks::DDPShocks, transfunc::Function,
# 				intdim::Type{T},
# 				tStateVectors::NTuple{2,Vector{Float64}}, tChoiceVectors::NTuple{1,Vector{Float64}},
# 				tOtherStateVectors::Tuple{Vector{Float64}},
# 				initializefunc::Function, initialize_exact::Bool) where T<:DDMIntDim
#
# 	nDim::Int64 = length(tStateVectors)
#
# 	# choose firms as last dimension because loop over firms (easier to code)
# 	dimShocks, nPeriods, nFirms = size(shocks.aSim)
# 	aSim = fill!(zeros(nDim, nPeriods+1, nFirms), NaN)
# 	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)
# 	mChoice = fill!(zeros(1, nPeriods+1, nFirms), NaN)
#
#
# 	#############
# 	# construct interpolator
# 	# if knots are not evenly spaced: need gridded interpolation
# 	itp_policy = interpolate(tStateVectors, sol.tmeshPolFun[1], Gridded(Linear()))
# 	itp_value = interpolate(tStateVectors, sol.meshValFun, Gridded(Linear()))
#
# 	# if knots are evenly spaced: scaled Bsplines (not working yet)
# 	# itp = interpolate(meshPolFun, BSpline(Linear()), OnGrid())
# 	# itp = scale(itp, vK, vAlog)
#
# 	# specify that if value requested outside of grid, return grid value
# 	# I am using my inbounds() function, but still there is an issue with values on the grid somehow..
# 	# hope to remove these lines in the future
# 	itp_policy = extrapolate(itp_policy, Flat())
# 	itp_value = extrapolate(itp_value, Flat())
#
# 	if initialize_exact
#         itp_policy0 = interpolate(tOtherStateVectors, sol.tmeshPolFunZero[1], Gridded(Linear()))
#         itp_policy0 = extrapolate(itp_policy0, Flat())
#     end
#
# 	# bounds such that choice variable does not exit bounds
# 	minChoice::Float64 = tChoiceVectors[1][1]
# 	maxChoice::Float64 = tChoiceVectors[1][end]
#
# 	# # for exogenous exit
#  	# do_exog_exit = isdefined(params,:π)
# 	# if do_exog_exit
# 	# 	do_exog_exit =  do_exog_exit * (params.π>0)
# 	# end
#
# 	# simulate one firm at a time, so that can do in parallel later
# 	# @sync @parallel for i = 1:nFirms
#
# 	for i = 1:nFirms
#
# 		vSim_i = fill!(zeros(nDim, nPeriods+1), NaN) # need NaN so that will register as exit
# 		vVal_i = fill!(zeros(nPeriods+1), NaN)
#
#         if initialize_exact
#             vSim_i[:,1] .= initializefunc(shocks.aInit[:,1,i], itp_policy0...)
#         else
#             vSim_i[:,1] .= initialize_simple(tStateVectors) # stable
#         end
#
# 		vVal_i[1] = itp_value(vSim_i[:,1]...)
#
# 		mShocks_i = shocks.aSim[:,:,i]
#
#
# 		for t = 1 : nPeriods
# 			# use interpolator to get optimal policy
# 			# policy = first state next period
# 			# should not be necessary if interpolator does it's job properly
#
# 			get_choice(vSim_it, itp_policy, minChoice, maxChoice) =
# 				inbounds(itp_policy(vSim_it...), minChoice, maxChoice)
#
# 			choice = get_choice(vSim_i[:,t], itp_policy, minChoice, maxChoice)
#
# 			# THIS BELOW LEADS TO TYPE INFERENCE PROBLEMS or is just slow
# 			# vSim_i[:,t+1] .= transfunc(SA, vSim_i[:,t], choice, mShocks_i[:,t])
#
# 			function fill_it!(vSim_it1, vSim_it, transfunc, intdim::Type{SA}, choice, mShocks_it)
# 				vSim_it1 .= transfunc(intdim, vSim_it, choice, mShocks_it)
# 			end
# 			function fill_it!(vSim_it1, vSim_it, transfunc, intdim::Type{intermediate}, choice, mShocks_it)
# 				vSim_it1[1] = choice
# 				vSim_it1[2:end] .= transfunc(intdim, vSim_it, mShocks_it)
# 			end
# 			function fill_it!(vSim_it1, vSim_it, transfunc, intdim::Type{separable}, choice, mShocks_it)
# 				vSim_it1[1] = choice
# 				vSim_it1[2:end] .= transfunc(intdim, vSim_it[2:end], mShocks_it)
# 			end
#
# 			fill_it!(vSim_i[:,t+1], vSim_i[:,t], transfunc, intdim, choice, mShocks_i[:,t])
#
#
# 			# if do_exog_exit
# 			# 	if i/nFirms > (1-params.π)^max(t-2, 0)
# 			# 		break
# 			# 	end
# 			# end # exit condition
#
# 		end
#
# 		aSim[:,:,i] = vSim_i
# 		mVal[1,:,i] = vVal_i
#
# 	end # firm loop
#
# 	return aSim
# 	# return cat(mVal, aSim, dims=1) #sdata(aSim) # convert shared array to standard array
# end # simulatemodel

# # for firm exit model
# function simulate(p::LearningKExit, sol::DDPSolution, shocks::DDPShocks)
# 	println("Note: still need to add computation of firm value in simulation of models with exit")
# 	# dmension of problem
# 	nDim = length(p.nNodes)
#
# 	# choose firms as last dimension because loop over firms (easier to code)
# 	 saSim = SharedArray{Float64}((nDim, nPeriods+1, nFirms))
# 	saSim[1,:,:] = NaN # need this to start iterator
# 	saSimExit = SharedArray{Float64}((1, nPeriods+1, nFirms))
#
#
# 	# interpolator for policy function
# 	itp = interpolate(p.tStateVectors, sol.meshPolFun, Gridded(Linear()))
#
# 	# interpolator for exit decision: will be 1 for exit, 0.5 for inbetween, otherwise zero
# 	itpExit = interpolate(p.tStateVectors, sol.meshExit, Gridded(Linear()))
#
# 	# interpolator for inital capital stock K0
# 	itp_K0 = interpolate(p.tStateVectors[2:3], sol.mChoice0, Gridded(Linear()))
#
# 	# interpolator for indirect value V0 (z, mu, P=P0)
# 	itp_V0 = interpolate(p.tStateVectors[1:3], sol.meshValFun[:,:,:,end], Gridded(Linear()))
#
# 	# V0 = (K = K0, z, mu, P = P0) indirect utility
# 	itp_V0ind = interpolate(p.tStateVectors[2:3], sol.mV0, Gridded(Linear()))
#
#
# 	# if do parallel, need to initialize a different Marsenne on each worker?
# 	RNG_Local = MersenneTwister(7) # create local RNG that does not interfere with BlackBoxOptim
# 	nShocks = size(p.mShocks,1)
#
# 	# simulate one firm at a time, so that can do in parallel later
# 	#@sync @parallel
# 	for f = 1:nFirms
#
# 		# iterate until a firm wants to enter
# 		while isnan(saSim[1,1,f])
# 			saSim[:, 1, f] = initialize(p, randn(RNG_Local, 3), itp_K0, itp_V0ind)
# 		end
#
# 		for t = 1 : nPeriods
#
# 			# policy = first state next period
# 			saSim[1, t+1, f] = itp(saSim[:,t,f]...)
#
# 			# other state vars
# 			saSim[2:end, t+1, f] = transfunc(p, saSim[2:end,t,f], shocks.aSim[:,t,f])
#
# 			# state vars predetermined: does firm want to exit next period?
# 			if itpExit[saSim[:, t+1, f]...] >= rand() # true if exit
# 				# if exit next period, we won't see its financial statement at end of year
# 				saSimExit[1,t+1,f] = 1.
#
# 				break
#
# 			elseif p.params.π >= rand()
# 				saSimExit[1,t+1,f] = 1.
#
# 				break
# 				# # already choose a new firm for then
# 				# # find new firm that enters
# 				# while true
# 				# 	saSim[:, t+1, f] = initialize(p, randn(RNG_Local, 3), itp_K0, itp_V0ind)
# 				# 	if !isnan(saSim[1, t+1, f]) # is available
# 				# 		break # break out of infinite loop
# 				# 	end
# 				# end
#
# 			end # exit condition
#
#
# 		end
#
#
# 	end # firm loop
#
# 	return vcat(sdata(saSim), sdata(saSimExit))  # convert shared array to standard array
# end # simulatemodel
