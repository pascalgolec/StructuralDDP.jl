
# function simulate(p::DDM)
# 	# seed!(42) # fix the seed
# 	simulate(p, drawshocks(p))
# end

# function simulate(p::DDM, sol::DDPSolution)
# 	# seed!(42) # fix the seed
# 	simulate(p, sol, drawshocks(p))
# end

simulate(p::DDM) = simulate(p, drawshocks(p))
simulate(p::DDM, sol::DDPSolution) = simulate(p, sol, drawshocks(p))
simulate(p::DDM, sol::DDPSolution, nPeriods::Int64, nFirms::Int64) =
	simulate(p, sol, drawshocks(p, nPeriods, nFirms))
simulate(p::DDM, shocks::DDPShocks) = simulate(p, solve(p), shocks)

function initialize_simple(p::DDM)
	# write a standard initialization function if want burn-in method
	# get the middle index of all state variables to initialize
	# convert(Vector{Float64},
	# 	getindex.(p.tStateVectors, Int.(floor.(length.(p.tStateVectors)./2))))
	getindex.(p.tStateVectors, Int.(floor.(length.(p.tStateVectors)./2)))
end

# perhaps better to return a tuple of arrays rather than an array
function simulate(p::DDM, sol::DDPSolution, shocks::DDPShocks;
        initialize_exact::Bool = sol.tmeshPolFunZero != nothing)

	# @unpack intdim, nPeriods, nFirms = p.params
	# transmethod = eval(intdim)
	# transmethod <: Union{SA,intermediate,separable} || error("intdim does not fit simulate")

	nDim = length(p.tStateVectors)
	nChoiceVars = length(sol.tmeshPolFun)

	# choose firms as last dimension because loop over firms (easier to code)
	dimShocks, nPeriods, nFirms = size(shocks.aSim)
	aSim = fill!(zeros(nDim, nPeriods+1, nFirms), NaN)
	mVal = fill!(zeros(1, nPeriods+1, nFirms), NaN)


	#############
	# construct interpolator
	# if knots are not evenly spaced: need gridded interpolation
	itp_policy = [interpolate(p.tStateVectors, polfun, Gridded(Linear())) for polfun in sol.tmeshPolFun]
	itp_value = interpolate(p.tStateVectors, sol.meshValFun, Gridded(Linear()))

	# if knots are evenly spaced: scaled Bsplines (not working yet)
	# itp = interpolate(meshPolFun, BSpline(Linear()), OnGrid())
	# itp = scale(itp, p.vK, p.vAlog)

	# specify that if value requested outside of grid, return grid value
	# I am using my inbounds() function, but still there is an issue with values on the grid somehow..
	# hope to remove these lines in the future
	itp_policy = [extrapolate(itp, Flat()) for itp in itp_policy]
	itp_value = extrapolate(itp_value, Flat())

	if initialize_exact
        # itp_policy0 = interpolate(p.tStateVectors[.!p.bEndogStateVars],
        #     sol.tmeshPolFunZero[1], Gridded(Linear()))
        itp_policy0 = [interpolate(p.tStateVectors[.!p.bEndogStateVars], polfun, Gridded(Linear())) for polfun in sol.tmeshPolFunZero]
        itp_policy0 = [extrapolate(itp, Flat()) for itp in itp_policy0]
    end


	# bounds such that choice variable does not exit bounds
	minChoice = [choicevec[1] for choicevec in p.tChoiceVectors]
	maxChoice = [choicevec[end] for choicevec in p.tChoiceVectors]

	# # for exogenous exit
 	# do_exog_exit = isdefined(p.params,:π)
	# if do_exog_exit
	# 	do_exog_exit =  do_exog_exit * (p.params.π>0)
	# end

	# simulate one firm at a time, so that can do in parallel later
	# @sync @parallel for i = 1:p.nFirms
	for i = 1:nFirms

		vSim_i = fill!(zeros(nDim, nPeriods+1), NaN) # need NaN so that will register as exit
		vVal_i = fill!(zeros(nPeriods+1), NaN)

        if initialize_exact
            vSim_i[:,1] .= initialize(p, shocks.aInit[:,1,i], itp_policy0...)
        else
            vSim_i[:,1] .= initialize_simple(p)
        end

		vVal_i[1] = itp_value(vSim_i[:,1]...)

		mShocks_i = shocks.aSim[:,:,i]

		for t = 1 : nPeriods
			# use interpolator to get optimal policy
			# policy = first state next period
			choice = [inbounds(itp_policy[i](vSim_i[:,t]...), minChoice[i], maxChoice[i])
				for i = 1:nChoiceVars]

			vSim_i[:,t+1] .= transfunc(p, SA, vSim_i[:,t], choice, mShocks_i[:,t])
			vVal_i[t+1] = itp_value(vSim_i[:,t+1]...)

			# if do_exog_exit
			# 	if i/nFirms > (1-p.params.π)^max(t-2, 0)
			# 		break
			# 	end
			# end # exit condition

		end

		aSim[:,:,i] = vSim_i
		mVal[1,:,i] = vVal_i

	end # firm loop

	return cat(mVal, aSim, dims=1) #sdata(aSim) # convert shared array to standard array
end # simulatemodel

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
