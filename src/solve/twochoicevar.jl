_solve(p::DDP{nStateVars,2}, method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
		disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
		rewardcall::Type{R}, monotonicity::Vector{Bool}, concavity::Vector{Bool}) where
			{nStateVars, T <: Separable_Union, R<:RewardCall} =
		_solve2(p.rewardfunc, method,
			mTransition, mReward,
			disp, disp_each_iter, max_iter, epsilon,
			rewardcall,
			monotonicity[1], monotonicity[2],
			concavity[1], concavity[2],
			p.tStateVectors,
			getchoicevars(p),
			getnonchoicevars(p),
			p.β)


function _solve2(rewardfunc, method::Type{T},
						mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
						disp::Bool, disp_each_iter::Int, max_iter::Int, epsilon::Float64,
						rewardcall::Type{R},
						monotonicity1::Bool, monotonicity2::Bool,
						concavity1::Bool, concavity2::Bool,
						tStateVectors,#::NTuple{2,Vector{Float64}},
						tChoiceVectors,
						tOtherStateVectors, #::NTuple{1,Vector{Float64}}
						β::Float64) where
						{T <: Separable_Union, R<:RewardCall}


    (nChoiceOne, nChoiceTwo) = length.(tChoiceVectors)
    (vChoiceOne, vChoiceTwo) = tChoiceVectors

	nChoices = nChoiceOne * nChoiceTwo
    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

	mValFun    = zeros((nChoices, nOtherStates)) # matrix so that can multiply with mTransition
	mValFunNew = zeros(nStates) # vector so that can index easily
	# mValFunNew = zeros((nChoices, nOtherStates)) # vector so that can index easily

    mPolFunInd1 = zeros(Int16, nStates)
    mPolFunInd2 = zeros(Int16, nStates)
    # mPolFunInd1 = zeros(Int16, nChoices, nOtherStates)
    # mPolFunInd2 = zeros(Int16, nChoices, nOtherStates)

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
	# mValFunDiff = zeros(nChoices, nOtherStates)
    # mValFunDiffAbs = zeros(nChoices, nOtherStates)

    liquidationvalue::Float64 = -Inf
    inactionvalue::Float64 = -Inf
    reward::Float64 = -Inf
	βEV::Float64 = -Inf

	valueHighSoFarOne::Float64 = -Inf
	valueProvisionalOne::Float64 = -Inf

	valueHighSoFarTwo::Float64 = -Inf
	valueProvisionalTwo::Float64 = -Inf

    iChoice1::Int16 = 0
    iChoice2::Int16 = 0

	# iChoice1Start::Int16 = 1
    # vChoice2Start = ones(Int64, nChoiceOne)
	# iChoice2inner = 0

	# iChoice2Start = 1
    vChoice2Start = ones(Int64, nChoiceOne)
    # vChoice1Start = ones(Int64, nChoiceTwo)
    vChoice1Start = 1
	iChoice1inner = 0

	i::Int64 = 0 # outer loop for other state vars
	cnt = initialize_counter()

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)

    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

        # @inbounds
		i = 0
		reset!(cnt)

		for ix in CartesianIndices(length.(tOtherStateVectors)) # other states
			i = i + 1
			cnt.i_exogstate += 1

            if monotonicity2
                vChoice2Start[:] .= 1
            end

            for l = 1:nChoiceTwo

				if monotonicity1
	                vChoice1Start = 1
	            end

                for j = 1:nChoiceOne # first state

                    valueHighSoFarTwo = -Inf
                    valueProvisionalTwo = -Inf
                    iChoice1 = 0
                    iChoice2 = 0

					cnt.i_state += 1
					cnt.i_choice = 0

					cnt.i_statechoice += vChoice2Start[j]-1 # compensate if start later
					cnt.i_choice += (vChoice2Start[j]-1)*nChoiceOne # compensate if start later

                    for lprime = vChoice2Start[j]:nChoiceTwo

        						# get optimal second choice variable conditional on first
        						valueHighSoFarOne = -Inf
        						valueProvisionalOne = -Inf
        						iChoice1inner  = 0

        						# find highest value for second state var
								cnt.i_choice += vChoice1Start-1
        						for jprime = vChoice1Start:nChoiceOne

									cnt.i_statechoice += 1
									cnt.i_choice += 1

            						# reward using pre_partial output matrix
            						if rewardcall == pre_partial
            							  reward = rewardfunc(mReward[cnt.i_state],
	  					                                     getindex.(tStateVectors, (j, l, ix.I...)),
														     (vChoiceOne[jprime], vChoiceTwo[lprime]))
            						elseif rewardcall == jit
            							  # need to be VERY careful with order of state vars here..
            							  reward = rewardfunc(getindex.(tStateVectors, (j, l, ix.I...)),
            												  (vChoiceOne[jprime], vChoiceTwo[lprime]))
            						elseif rewardcall == pre
            							  # reward = mReward[jprime + nChoiceOne * (lprime-1), cnt.i_state] # nChoices x nStates
            							  reward = mReward[cnt.i_choice, cnt.i_state] # nChoices x nStates
            						end

					                if method == Separable_ExogStates # mβEV is nChoices x nStochStates
											# βEV = mβEV[jprime + nChoiceOne * (lprime-1), cnt.i_exogstate]
											βEV = mβEV[cnt.i_choice, cnt.i_exogstate]
        							elseif method == Separable_States # mβEV is nChoices x nStates
            								# βEV = mβEV[jprime + nChoiceOne * (lprime-1), cnt.i_state]
            								βEV = mβEV[cnt.i_choice, cnt.i_state]
        							else method == Separable # mβEV is nChoices x (nChoices * nStates)
            								# βEV = mβEV[jprime + nChoiceOne * (lprime-1), cnt.i_statechoice]
            								βEV = mβEV[cnt.i_choice, cnt.i_statechoice]
        							end

									# if i==1
									# 	@show cnt.i_choice
									# 	@show jprime + nChoiceOne * (lprime-1)
									# end
									# if cnt.i_choice != jprime + nChoiceOne * (lprime-1)
									# 	@show jprime
									# 	@show lprime
									# 	@show nChoiceOne
									# 	@show vChoice2Start[j]
									# 	@show cnt.i_choice
									# 	@show jprime + nChoiceOne * (lprime-1)
									# 	error("stop")
									# end

        							valueProvisionalOne = reward + βEV

        							if valueProvisionalOne >= valueHighSoFarOne
            							   valueHighSoFarOne = valueProvisionalOne
            							   iChoice1inner = jprime
			                        elseif concavity1
            							break
        							end

        						end # jprime

                        valueProvisionalTwo = valueHighSoFarOne

                        if (valueProvisionalTwo>=valueHighSoFarTwo)
                            valueHighSoFarTwo = valueProvisionalTwo
                            iChoice1 = iChoice1inner
                            iChoice2 = lprime
                        elseif concavity2
                            break # break when we have achieved the max
                        end

                    end #lprime

                    # if isdefined(p.params, :F)
                    #     (inactionvalue, iChoice2inaction) =
                    #                     getsecondchoice(p,
                    #                                     passive, # whether capital choice active or passive
                    #                                     mReward, mβEV, i, j, l,
                    #                                     nChoiceOne, vChoiceOne, j, # inaction: jprime = j
                    #                                     nChoiceTwo, vChoiceTwo,
                    #                                     vChoice2Start[j])
                    #
                    #     if valueHighSoFarOne <= inactionvalue
                    #         # don't have interior K'
                    #         valueHighSoFarOne = inactionvalue
                    #         iChoice1 = j
                    #         iChoice2 = iChoice2inaction
                    #     end
                    # end

                    mValFunNew[cnt.i_state] = valueHighSoFarTwo
                    mPolFunInd1[cnt.i_state] = iChoice1
                    mPolFunInd2[cnt.i_state] = iChoice2
                    # mValFunNew[j + nChoiceOne * (l-1), i] = valueHighSoFarTwo
                    # mPolFunInd1[j + nChoiceOne * (l-1), i] = iChoice1
                    # mPolFunInd2[j + nChoiceOne * (l-1), i] = iChoice2

					if monotonicity1
                        vChoice1Start = iChoice1
                    end
					if monotonicity2
					  vChoice2Start[j] = iChoice2
					end

                end # j (first choice)

            end # l (second choice)

        end #i (stochastic variables)

		# @show size(mValFunNew)
		# @show size(mValFun)
		# @show size(mValFun[:])
        # mValFunDiff .= mValFunNew .- mValFun
        mValFunDiff .= mValFunNew .- @view(mValFun[:])
        mValFunDiffAbs .= abs.(mValFunDiff)
        maxDifference  = maximum(mValFunDiffAbs)

        # mcqueen!(mValFunNew, mValFun, β) #--> does not even converge somehow

        copyto!(mValFun, mValFunNew)

        # howards improvement?

        if disp
            display_iter(iteration, disp_each_iter, maxDifference)
        end

		iteration += 1

        if iteration > max_iter
            println("WARNING: maximum iterations exceeded")
            # println("Parameters used:")
			# println(p.params)
            break
        end

    end #while


    # println("Iteration ", string(iteration), ": Sup change is = ", string(maxDifference))

    mPolFun1 = vChoiceOne[mPolFunInd1]
    mPolFun2 = vChoiceTwo[mPolFunInd2]

    nNodes = length.(tStateVectors)
    meshPolFun1 = reshape(mPolFun1, tuple(nNodes...))
    meshPolFun2 = reshape(mPolFun2, tuple(nNodes...))

    meshValFun = reshape(mValFun, tuple(nNodes...))
    # meshExit   = reshape(mExit, tuple(nNodes...))

	return meshValFun, (meshPolFun1, meshPolFun2)
end # solve
