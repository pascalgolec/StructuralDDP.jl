_solve(p::DDP{nStateVars,2}, method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
		disp::Bool, rewardmat::Symbol, monotonicity::Vector{Bool}, concavity::Vector{Bool}) where
			{nStateVars, T <: Separable_Union} =
		_solve2(p.rewardfunc, method,
			mTransition, mReward,
			disp, rewardmat,
			monotonicity[1], monotonicity[2],
			concavity[1], concavity[2],
			p.tStateVectors,
			getchoicevars(p),
			getnonchoicevars(p),
			p.β)


function _solve2(rewardfunc::Function, method::Type{T},
						mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
						disp::Bool, rewardmat,
						monotonicity1::Bool, monotonicity2::Bool,
						concavity1::Bool, concavity2::Bool,
						tStateVectors,#::NTuple{2,Vector{Float64}},
						tChoiceVectors,
						tOtherStateVectors, #::NTuple{1,Vector{Float64}}
						β::Float64) where
						T <: Separable_Union


    (nChoiceOne, nChoiceTwo) = length.(tChoiceVectors)
    (vChoiceOne, vChoiceTwo) = tChoiceVectors

	nChoices = nChoiceOne * nChoiceTwo
    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

    mValFun    = zeros((nChoices, nOtherStates))
    mValFunNew = zeros((nChoices, nOtherStates))

    mPolFunInd1 = zeros(Int16, nChoices, nOtherStates)
    mPolFunInd2 = zeros(Int16, nChoices, nOtherStates)

    mβEV = zeros(nChoices, size(mTransition, 1)) # depends on intdim

    # VFI initialization
    maxDifference::Float64 = 10000.0 # to initialize
    tolerance = 1.E-8
    iteration::Int64 = 0 # initialize counter

    # initialize for less memory allocation
	mValFunDiff = zeros(nChoices, nOtherStates)
    mValFunDiffAbs = zeros(nChoices, nOtherStates)

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

	iChoice1Start::Int16 = 1
    vChoice2Start = ones(Int64, nChoiceOne)

	i::Int64 = 0 # outer loop for other state vars

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)


    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

        # @inbounds
		i = 0
		for ix in CartesianIndices(length.(tOtherStateVectors)) # other states
			i = i + 1

            if monotonicity2
                vChoice2Start[:] .= 1
            end

            for l = 1:nChoiceTwo

                # We start from previous choice (monotonicity of policy function)
                if monotonicity1
                   iChoice1Start = 1
                end

                for j = 1:nChoiceOne # first state

                    valueHighSoFarOne = -Inf
                    valueProvisionalOne = -Inf
                    iChoice1 = 0
                    iChoice2 = 0

                    for jprime = iChoice1Start:nChoiceOne

        						# get optimal second choice variable conditional on first
        						valueHighSoFarTwo = -Inf
        						valueProvisionalTwo = -Inf
        						iChoice2inner  = 0

        						# find highest value for second state var
        						for lprime = vChoice2Start[j]:nChoiceTwo

            						# reward using prebuild_partial output matrix
            						if rewardmat == :prebuild_partial
            							  reward = rewardfunc(mReward[j + nChoiceOne*(l-1),i],
            												  vChoiceOne[j], vChoiceOne[jprime],
            												  vChoiceTwo[l], vChoiceTwo[lprime])
            						elseif rewardmat == :nobuild
            							  # need to be VERY careful with order of state vars here.. could get fucked up..
            							  reward = rewardfunc(getindex.(tStateVectors, [j, l, ix.I...]),
            												  [vChoiceOne[jprime], vChoiceTwo[lprime]])
            						elseif rewardmat == :prebuild
            							  # jprime is first choice var, changes faster
            							  reward = mReward[jprime + nChoiceOne * (lprime-1), j + nChoiceOne * (l-1) + (nChoiceOne*nChoiceTwo) * (i-1)] # nChoices x nStates
            						end

					                if method == Separable_ExogStates # mβEV is nChoices x nStochStates
            								βEV = mβEV[jprime + nChoiceOne * (lprime-1), i] # mβEV is already discounted, jprime = j for inactive
        							elseif method == Separable_States # mβEV is nChoices x nStates
            								βEV = mβEV[jprime + nChoiceOne * (lprime-1), j + nChoiceOne * (l-1) + (nChoiceOne*nChoiceTwo) * (i-1)]
        							else method == Separable # mβEV is mβEV is nChoices x (nChoices * nStates)
            								# error("separable not done yet")
            								βEV = mβEV[jprime + nChoiceOne * (lprime-1),
            								  jprime + nChoiceOne * (lprime-1) + (nChoiceOne*nChoiceTwo) *
            									  (-1 + j + nChoiceOne * (l-1) + (nChoiceOne*nChoiceTwo) * (i-1)) ]
        							end

        							valueProvisionalTwo = reward + βEV

        							if valueProvisionalTwo >= valueHighSoFarTwo
            							   valueHighSoFarTwo = valueProvisionalTwo
            							   iChoice2inner = lprime
			                        elseif concavity2
            							break
        							end

        						end # lprime

        						# return valueHighSoFarTwo, iChoice2

                        valueProvisionalOne = valueHighSoFarTwo

                        if (valueProvisionalOne>=valueHighSoFarOne)
                            valueHighSoFarOne = valueProvisionalOne
                            iChoice1 = jprime
                            iChoice2 = iChoice2inner
                            if monotonicity1
                	          iChoice1Start = jprime
                            end
                        elseif concavity1
                            break # We break when we have achieved the max
                        end

                    end #jprime

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

                    mValFunNew[j + nChoiceOne * (l-1), i] = valueHighSoFarOne
                    mPolFunInd1[j + nChoiceOne * (l-1), i] = iChoice1
                    mPolFunInd2[j + nChoiceOne * (l-1), i] = iChoice2
                    if monotonicity2
                        vChoice2Start[j] = iChoice2
                    end

                end #j (capital)

            end # l (second choice)

        end #i (stochastic variables)

        mValFunDiff .= mValFunNew .- mValFun
        mValFunDiffAbs .= abs.(mValFunDiff)
        maxDifference  = maximum(mValFunDiffAbs)

        # mcqueen!(mValFunNew, mValFun, β) #--> does not even converge somehow

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

    mPolFun1 = vChoiceOne[mPolFunInd1]
    mPolFun2 = vChoiceTwo[mPolFunInd2]

    nNodes = length.(tStateVectors)
    meshPolFun1 = reshape(mPolFun1, tuple(nNodes...))
    meshPolFun2 = reshape(mPolFun2, tuple(nNodes...))

    meshValFun = reshape(mValFun, tuple(nNodes...))
    # meshExit   = reshape(mExit, tuple(nNodes...))

	return meshValFun, (meshPolFun1, meshPolFun2)
end # solve
