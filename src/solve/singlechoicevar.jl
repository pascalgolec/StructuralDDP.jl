########################
##### Solver ######
########################

_solve(p::DiscreteDynamicProblem{nStateVars, 1, T},
		method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
		disp::Bool, rewardmat::Symbol, monotonicity::Bool, concavity::Bool) where
			{T <: Union{separable, intermediate}, nStateVars} =
		_solve1(p.rewardfunc, method,
			mTransition, mReward,
			disp, rewardmat, monotonicity, concavity,
			p.tStateVectors, p.tStateVectors[p.bEndogStateVars][1], p.tStateVectors[.!p.bEndogStateVars],
			# p.bEndogStateVars,
			p.params.β)

function _solve1(rewardfunc::Function, method::Type{T},
				mTransition::Array{Float64,2}, mReward::Array{Float64,2},
				disp::Bool, rewardmat::Symbol, monotonicity::Bool, concavity::Bool,
				tStateVectors,#::NTuple{2,Vector{Float64}},
				vChoices::Vector{Float64},
				tOtherStateVectors, #::NTuple{1,Vector{Float64}}
				β::Float64) where
				T <: Union{separable, intermediate}

    nChoices = length(vChoices)

    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

    mValFun    = zeros((nChoices, nOtherStates))
	mValFunNew = zeros((nChoices, nOtherStates))

	mPolFunInd = zeros(Int16, nChoices, nOtherStates)

    mβEV = zeros(nChoices, size(mTransition, 1)) # depends on intdim

    # VFI initialization
    maxDifference::Float64 = 10000.0 # to initialize
    tolerance = 1.E-8
    iteration::Int64 = 0 # initialize counter

    # initialize for less memory allocation
	mValFunDiff = zeros(nChoices,nOtherStates)
    mValFunDiffAbs = zeros(nChoices,nOtherStates)

    liquidationvalue::Float64 = -10000.
    inactionvalue::Float64 = -10000.
    reward::Float64 = -10000.

	valueHighSoFar::Float64 = -1000.
	valueProvisional::Float64 = -1000.
    iChoiceStart::Int16 = 1
    iChoice::Int16 = 1
	i::Int64 = 0 # outer loop for other state vars

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)

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

					# if typeof(p) == NeoClassicalVolatilityAge
					# 	age = mgridstateother[i,2]
                    # 	reward = rewardfunc(p, mReward[j,i], vChoices[j], vChoices[jprime], age)
					# elseif typeof(p) == LearningKVolatilityAge
					# 	age = mgridstateother[i,3]
                    # 	reward = rewardfunc(p, mReward[j,i], vChoices[j], vChoices[jprime], age)
					# elseif typeof(p) == NewIdeasAge
					# 	age = mgridstateother[i,3]
                    # 	reward = rewardfunc(p, mReward[j,i], vChoices[j], vChoices[jprime], age)
					# else
                    # 	reward = rewardfunc(p, mReward[j,i], vChoices[j], vChoices[jprime])
					# end

					# reward using prebuild_partial output matrix
					if rewardmat == :prebuild_partial
						reward = rewardfunc(mReward[j,i], vChoices[j], vChoices[jprime])
					elseif rewardmat == :nobuild
						# need to be VERY careful with order of state vars here.. could get fucked up..
						# reward = rewardfunc(p, getindex.(tStateVectors, [j, ix.I...]), vChoices[jprime])
						reward = rewardfunc(getindex.(tStateVectors, [j, ix.I...]), vChoices[jprime])
					elseif rewardmat == :prebuild
						reward = mReward[jprime, j + nChoices * (i-1)] # nChoices x nStates
					end

                    if method == separable
	                    valueProvisional = reward + mβEV[jprime, i] # mβEV is already discounted
                    else # method == intermediate
                        valueProvisional = reward + mβEV[jprime, j + nChoices *(i-1)] # mβEV is already discounted
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

				# # compare with discontinuity
                # if isdefined(p.params, :F) || typeof(p) <: SalesAdjCosts
                #     # compare with inaction: investing enough to keep K constant
    			# 	reward = rewardfuncinaction(p, mReward[j,i], vChoices[j])
				#
                #     if intdim == :separable
                #         inactionvalue = reward + mβEV[j, i]
                #     elseif intdim == :intermediate
                #         inactionvalue = reward + mβEV[j, j + nChoices *(i-1)] # mβEV is already discounted
                #     end
				#
                #     if valueHighSoFar <= inactionvalue
                #         # don't have interior K'
                #         valueHighSoFar = inactionvalue
                #         iChoice = j
                #     end
                # end


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
    	    if mod(iteration,10)==0 || iteration == 1
    	        println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
    	    end
        end

        if iteration > 1000
            println("WARNING: maximum iterations exceeded")
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
