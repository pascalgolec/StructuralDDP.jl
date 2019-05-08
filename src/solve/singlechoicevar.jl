mutable struct Counters
	"""linear counter for Separable_ExogStates"""
	i_exogstate::Int64
	"""linear counter for Separable_States"""
	i_state::Int64
	"""linear counter for Separable"""
	i_statechoice::Int64
	"""linear counter for choice"""
	i_choice::Int64
end

initialize_counter() = Counters(0,0,0,0)
getcounter(cnt::Counters, intdim::Type{Separable_ExogStates}) = cnt.i_exogstate
getcounter(cnt::Counters, intdim::Type{Separable_States}) = cnt.i_state
getcounter(cnt::Counters, intdim::Type{Separable}) = cnt.i_statechoice
function reset!(cnt::Counters)
	cnt.i_exogstate = 0
	cnt.i_state = 0
	cnt.i_statechoice = 0
	cnt.i_choice = 0
end

getreward(rewardcall::Type{pre_partial}, rewardfunc, mReward, tStateVectors, vChoices,
	cnt, ix, j, jprime) =
	rewardfunc(mReward[cnt.i_state], getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])

getreward(rewardcall::Type{jit}, rewardfunc, mReward, tStateVectors, vChoices,
	cnt, ix, j, jprime) =
	rewardfunc(getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])

getreward(rewardcall::Type{pre}, rewardfunc, mReward, tStateVectors, vChoices,
	cnt, ix, j, jprime) =  mReward[jprime, cnt.i_state] # nChoices x nStates


# function getreward(rewardcall, rewardfunc, mReward, tStateVectors, vChoices,
# 	cnt, ix, j, jprime)
# 	if rewardcall == :pre_partial
# 		return rewardfunc(mReward[cnt.i_state], getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])
# 	elseif rewardcall == :jit
# 		# need to be VERY careful with order of state vars here..
# 		return rewardfunc(getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])
# 	elseif rewardcall == :pre
# 		return mReward[jprime, cnt.i_state] # nChoices x nStates
# 	end
# end

########################
##### Solver ######
########################
_solve(p::DDP{NS,1},
	mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
	opts::SolverOptions) where NS =
	_solve(p, mTransition, mReward, opts, p.tStateVectors,
		getchoicevars(p.tStateVectors, p.tChoiceVectors)[1],
		getnonchoicevars(p.tStateVectors, p.tChoiceVectors),)

# function barrier for tuples
function _solve(p::DDP{NS,1},
				mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
				opts::SolverOptions,
				tStateVectors,
				vChoices::Vector{Float64},
				tOtherStateVectors) where NS

	@unpack β, rewardfunc = p
	@unpack disp, disp_each_iter, max_iter, epsilon, rewardcall,
		monotonicity, concavity, intdim = opts
	get_additional_index = p.options.get_additional_index

    nChoices = length(vChoices)

    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

    mValFun    = zeros((nChoices, nOtherStates)) # matrix so that can multiply with mTransition
	mValFunNew = zeros(nStates) # vector so that can index easily

	mPolFunInd = zeros(Int16, nStates)

    mβEV = zeros(nChoices, size(mTransition, 1)) # depends on intdim

    # VFI initialization
    maxDifference::Float64 = Inf # to initialize
    iteration::Int64 = 0 # initialize counter
	if β == 0.0
        tolerance = Inf
    elseif β == 1.0
        throw(ArgumentError("method invalid for beta = 1"))
    else
        tolerance = epsilon * (1-β) / (2*β)
    end

    # initialize for less memory allocation
	mValFunDiff = zeros(nStates)
    mValFunDiffAbs = zeros(nStates)

    # liquidationvalue::Float64 = -10000.
    # inactionvalue::Float64 = -10000.
    reward::Float64 = -Inf

	valueHighSoFar::Float64 = -Inf
	valueProvisional::Float64 = -Inf
    iChoiceStart::Int16 = 1
    iChoice::Int16 = 1

	cnt = initialize_counter()

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)

	# whether to check outside of monotonicity/concavity
	try_additional = any((monotonicity, concavity)) && get_additional_index!=nothing

	check_jprime  = -1

    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

		reset!(cnt)

		# @inbounds
		for ix in CartesianIndices(length.(tOtherStateVectors)) # other states

			cnt.i_exogstate += 1

            if monotonicity
               iChoiceStart = 1
            end

	        for j = 1:nChoices # first state

	            valueHighSoFar = -Inf
	            iChoice  = 1

				cnt.i_state += 1

				cnt.i_statechoice += iChoiceStart-1 # compensate if start later

	            for jprime = iChoiceStart:nChoices

					cnt.i_statechoice += 1

					reward = getreward(rewardcall, rewardfunc, mReward, tStateVectors, vChoices,
						cnt, ix, j, jprime)

                    valueProvisional = reward + mβEV[jprime, getcounter(cnt, intdim)]

	                if (valueProvisional>=valueHighSoFar)
    	            	valueHighSoFar = valueProvisional
    	            	iChoice = jprime
	                elseif concavity
						cnt.i_statechoice += nChoices - jprime # adjust counter if finish early
	            	    break # break bc value will only be lower from now on
					end

	            end #jprime

				if monotonicity
				  iChoiceStart = iChoice
				end

				if try_additional
					# check one more value outside of monotonicity and conc
					check_jprime = get_additional_index((j, ix.I...))

					cnt.i_statechoice += check_jprime - nChoices

					reward = getreward(rewardcall, rewardfunc, mReward, tStateVectors, vChoices,
						cnt, ix, j, check_jprime)

					valueProvisional = reward + mβEV[check_jprime, getcounter(cnt, intdim)]

					if valueHighSoFar <= valueProvisional # check is better than interior
						valueHighSoFar = valueProvisional
						iChoice = check_jprime
					end

					cnt.i_statechoice += nChoices - check_jprime

				end

                # # compare with liquidation
                # if typeof(p) == LearningKExit
                #     liquidationvalue = mReward[j,i] + p.β*p.κ*(1-p.δ)*vChoices[j]
                #
				# 	# liquidation value is when desinvest to zero, don't need to pay future fixed costs
				# 	# can not be smaller than zero because of limited liability
				# 	# liquidationvalue = max(0, rewardfunc(mReward[j,i], vChoices[j], 0., p))
                #
				# 	# liquidationvalue = dividends(rewardfunc(mReward[j,i], vChoices[j], 0., p), p)
                #     if valueHighSoFar < liquidationvalue
                #        valueHighSoFar = liquidationvalue
                #        mExit[j,i] = true
                #        #iChoice = j dont need this, its just for iterating
                #     else
                #        mExit[j,i] = false
                #     end
                # end

	            mValFunNew[cnt.i_state] = valueHighSoFar
	            mPolFunInd[cnt.i_state] = iChoice

	        end #j

	    end #i

        mValFunDiff .= mValFunNew .- @view(mValFun[:])
        mValFunDiffAbs .= abs.(mValFunDiff)
	    maxDifference  = maximum(mValFunDiffAbs)

        # mcqueen!(mValFunNew, mValFun, p.β)

        copyto!(mValFun, mValFunNew)

        # howards improvement?

        if disp
			display_iter(iteration, disp_each_iter, maxDifference)
        end

		iteration += 1

        if iteration > max_iter
            println("WARNING: $max_iter maximum iterations exceeded")
            # println("Parameters used:")
			# println(p.params)
            break
        end



	end #while


    # println("Iteration ", string(iteration), ": Sup change is = ", string(maxDifference))

    mPolFun = vChoices[mPolFunInd]

    nNodes = length.(tStateVectors)
    meshPolFun = reshape(mPolFun, tuple(nNodes...))
    meshValFun = reshape(mValFun, tuple(nNodes...))
    # meshExit   = reshape(mExit, tuple(p.nNodes...))

    meshValFun, (meshPolFun,)

end # solve

# function mcqueen!(v, vold, beta)
#     b_l = beta /(1-beta) * minimum(v-vold);
#     b_u = beta /(1-beta) * maximum(v-vold);
#     # @show (b_l+b_u)/2
#     v[:] = v[:] + (b_l+b_u)/2;
#     nothing
# end
