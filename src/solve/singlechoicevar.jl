########################
##### Solver ######
########################
_solve(p::DDP{nStateVars,1},
		method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
		disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
		rewardcall::Symbol, monotonicity::Bool, concavity::Bool) where
			{nStateVars, T <: Separable_Union} =
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
				mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
				disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
				rewardcall::Symbol, monotonicity::Bool, concavity::Bool,
				tStateVectors,#::NTuple{2,Vector{Float64}},
				vChoices::Vector{Float64},
				tOtherStateVectors, #::NTuple{1,Vector{Float64}}
				β::Float64,
				get_additional_index) where
				T <: Separable_Union

    nChoices = length(vChoices)

    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

    mValFun    = zeros((nChoices, nOtherStates))
	mValFunNew = zeros((nChoices, nOtherStates))

	mPolFunInd = zeros(Int16, nChoices, nOtherStates)

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
	mValFunDiff = zeros(nChoices,nOtherStates)
    mValFunDiffAbs = zeros(nChoices,nOtherStates)

    # liquidationvalue::Float64 = -10000.
    # inactionvalue::Float64 = -10000.
    reward::Float64 = -Inf

	valueHighSoFar::Float64 = -Inf
	valueProvisional::Float64 = -Inf
    iChoiceStart::Int16 = 1
    iChoice::Int16 = 1
	i::Int64 = 0 # outer loop for other state vars

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
		# for i = 1:nOtherStates # other states
		i = 0
		for ix in CartesianIndices(length.(tOtherStateVectors)) # other states
			i = i + 1

	        # We start from previous choice (monotonicity of policy function)
            if monotonicity
               iChoiceStart = 1
            end

	        for j = 1:nChoices # first state

	            valueHighSoFar = -1000.0
	            iChoice  = 1

	            for jprime = iChoiceStart:nChoices

					if rewardcall == :pre_partial
						reward = rewardfunc(mReward[j,i], getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])
					elseif rewardcall == :jit
						# need to be VERY careful with order of state vars here..
						reward = rewardfunc(getindex.(tStateVectors, (j, ix.I...)), vChoices[jprime])
					elseif rewardcall == :pre
						reward = mReward[jprime, j + nChoices * (i-1)] # nChoices x nStates
					end

                    if method == Separable_ExogStates
	                    valueProvisional = reward + mβEV[jprime, i] # mβEV is nChoices x nExogStates
                    elseif method == Separable_States
						valueProvisional = reward + mβEV[jprime, j + nChoices *(i-1)] # mβEV is nChoices x nStates
                    else # method == Separable
						valueProvisional = reward + mβEV[jprime, jprime + nChoices*(j-1 + nChoices *(i-1))] # mβEV is nChoices x (nStates * nChoices)
                    end

	                if (valueProvisional>=valueHighSoFar)
    	            	valueHighSoFar = valueProvisional
    	            	iChoice = jprime
                        if monotonicity
            	          iChoiceStart = jprime
                        end
	                elseif concavity
	            	    break # We break when we have achieved the max
					end

	            end #jprime

				if try_additional

					# check one more value outside of monotonicity and conc
					check_jprime = get_additional_index((j, ix.I...))

					if rewardcall == :pre_partial
						reward = rewardfunc(mReward[j,i], getindex.(tStateVectors, (j, ix.I...)), vChoices[check_jprime])
					elseif rewardcall == :jit
						# need to be VERY careful with order of state vars here..
						reward = rewardfunc(getindex.(tStateVectors, (j, ix.I...)), vChoices[check_jprime])
					elseif rewardcall == :pre
						reward = mReward[check_jprime, j + nChoices * (i-1)] # nChoices x nStates
					end

					if method == Separable_ExogStates
	                    valueProvisional = reward + mβEV[check_jprime, i] # mβEV is nChoices x nExogStates
                    elseif method == Separable_States
						valueProvisional = reward + mβEV[check_jprime, j + nChoices *(i-1)] # mβEV is nChoices x nStates
                    else # method == Separable
						valueProvisional = reward + mβEV[check_jprime, check_jprime + nChoices*(j-1 + nChoices *(i-1))] # mβEV is nChoices x (nStates * nChoices)
                    end

					if valueHighSoFar <= valueProvisional # check is better than interior
						valueHighSoFar = valueProvisional
						iChoice = check_jprime
					end

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

	            mValFunNew[j,i] = valueHighSoFar
	            mPolFunInd[j,i] = iChoice

	        end #j

	    end #i

        mValFunDiff .= mValFunNew .- mValFun
        mValFunDiffAbs .= abs.(mValFunDiff)
	    maxDifference  = maximum(mValFunDiffAbs)

        # mcqueen!(mValFunNew, mValFun, p.β) #--> does not even converge somehow

        copyto!(mValFun, mValFunNew)

        # howards improvement?

	    iteration += 1

        if disp
			display_iter(iteration, disp_each_iter, maxDifference)
        end

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
