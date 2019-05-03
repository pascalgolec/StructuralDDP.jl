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
getcounter(cnt::Counters, method::Type{Separable_ExogStates}) = cnt.i_exogstate
getcounter(cnt::Counters, method::Type{Separable_States}) = cnt.i_state
getcounter(cnt::Counters, method::Type{Separable}) = cnt.i_statechoice
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
_solve(p::DDP{nStateVars,1},
		method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
		disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
		rewardcall::Type{R}, monotonicity::Bool, concavity::Bool) where
			{nStateVars, T <: Separable_Union, R<:RewardCall} =
		_solve1(p.rewardfunc, method,
			mTransition, mReward,
			disp, disp_each_iter, max_iter, epsilon,
			rewardcall, monotonicity, concavity,
			p.tStateVectors,
			getchoicevars(p.tStateVectors, p.tChoiceVectors)[1],
			getnonchoicevars(p.tStateVectors, p.tChoiceVectors),
			p.β,
			p.options.get_additional_index)


# precondition is separable. have EndogStatevectors and exogstatevectors
function _solve1(rewardfunc, method::Type{T},
				mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
				disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
				rewardcall::Type{R}, monotonicity::Bool, concavity::Bool,
				tStateVectors,#::NTuple{2,Vector{Float64}},
				vChoices::Vector{Float64},
				tOtherStateVectors, #::NTuple{1,Vector{Float64}}
				β::Float64,
				get_additional_index) where
				{T <: Separable_Union, R<:RewardCall}

    nChoices = length(vChoices)

    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

    mValFun    = zeros((nChoices, nOtherStates)) # matrix so that can multiply with mTransition
	mValFunNew = zeros(nStates) # vector so that can index easily
	# mValFunNew = zeros((nChoices, nOtherStates)) # vector so that can index easily

	mPolFunInd = zeros(Int16, nStates)
	# mPolFunInd = zeros(Int16, nChoices, nOtherStates)

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
	# mValFunDiff = zeros(nChoices,nOtherStates)
    # mValFunDiffAbs = zeros(nChoices,nOtherStates)

    # liquidationvalue::Float64 = -10000.
    # inactionvalue::Float64 = -10000.
    reward::Float64 = -Inf

	valueHighSoFar::Float64 = -Inf
	valueProvisional::Float64 = -Inf
    iChoiceStart::Int16 = 1
    iChoice::Int16 = 1


	# need different counters depending on integration dimension and rewardcall
	# i::Int64 = 0 # linear counter for exogenous state var
	# i_state::Int64 = 0 # linear counter for state var
	# i_statechoice::Int64 = 0 # linear counter for state*choice var
	cnt = initialize_counter()

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)

	# whether to check outside of monotonicity/concavity
	try_additional = any((monotonicity, concavity)) && get_additional_index!=nothing

	check_jprime  = -1


    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

        # inbounds does help by 50%
	    # @inbounds
		reset!(cnt)

		for ix in CartesianIndices(length.(tOtherStateVectors)) # other states

			cnt.i_exogstate += 1

	        # We start from previous choice (monotonicity of policy function)
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

					# somehow it's slow when I put it into function.. should have dispatch for rewardcall
					reward = getreward(rewardcall, rewardfunc, mReward, tStateVectors, vChoices,
						cnt, ix, j, jprime)

                    valueProvisional = reward + mβEV[jprime, getcounter(cnt, method)]

	                if (valueProvisional>=valueHighSoFar)
    	            	valueHighSoFar = valueProvisional
    	            	iChoice = jprime
                        if monotonicity
            	          iChoiceStart = jprime
                        end
	                elseif concavity
						cnt.i_statechoice += nChoices - jprime # adjust counter if finish early
	            	    break # We break when we have achieved the max
					end

	            end #jprime

				if try_additional
					# check one more value outside of monotonicity and conc
					check_jprime = get_additional_index((j, ix.I...))

					cnt.i_statechoice += check_jprime - nChoices

					reward = getreward(rewardcall, rewardfunc, mReward, tStateVectors, vChoices,
						cnt, ix, j, check_jprime)

					valueProvisional = reward + mβEV[check_jprime, getcounter(cnt, method)]

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
	            # mValFunNew[j,cnt.i_exogstate] = valueHighSoFar
	            mPolFunInd[cnt.i_state] = iChoice
	            # mPolFunInd[j,cnt.i_exogstate] = iChoice

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
