
function getsecondchoice(p::TwoChoiceVar,
                         # action::Type{act},
                         mReward::Union{Array{Float64,2},Nothing},
                         mβEV::Array{Float64,2},
                         i::Int64, j::Int64, l::Int64,
                         nChoiceOne::Int64, vChoiceOne::Vector{Float64}, jprime::Int64,
                         nChoiceTwo::Int64, vChoiceTwo::Vector{Float64},
                         iChoice2Start::Int64,
						 ixI::NTuple{N,Int64} where N,
						 intdim::Type{T},
						 rewardmat::Symbol,
						 concavity::Vector{Bool}) where
						 	T <: Union{separable, intermediate}
                            # {act<:capitalaction}


    # get optimal second choice variable conditional on first
    valueHighSoFarTwo = -Inf
    valueProvisionalTwo = -Inf
    iChoice2  = 0

    # iChoice2Start = vChoice2Start[j]

    # find highest value for second state var
    for lprime = iChoice2Start:nChoiceTwo

        # if typeof(p) <: FinancingDebt &&
        #         (1+p.r*(1-p.τ))*vChoiceTwo[lprime]>p.θ*(1-p.δ)*vChoiceOne[jprime]
        #         # collateral constraint violated
        #         valueProvisionalTwo = -Inf
        #         break
        # end

        # if typeof(p) <: Intangible
        #     output = mReward[j + nChoiceOne*(l-1), i] # first dimension of mReward is nChoiceOne x nChoiceTwo if second variable affects output
        # elseif typeof(p) <: Union{IEI, IEIR, FinancingLeverageSimple}
        #     output = mReward[j,i] # second state variable does not affect output: output matrix is smaller
        # end

        # if action == active
		#
        #     reward = rewardfunc(p, output, vChoiceOne[j], vChoiceOne[jprime],
        #                         vChoiceTwo[l], vChoiceTwo[lprime])
		#
        # elseif action == passive
        #     # second state variable does not affect output
        #     reward = rewardfuncinaction(p, output, vChoiceOne[j],
        #                                 vChoiceTwo[l], vChoiceTwo[lprime])
        # end

		# preallocate reward?

		# reward using prebuild_partial output matrix
		if rewardmat == :prebuild_partial
			reward = rewardfunc(p, mReward[j + nChoiceOne*(l-1),i],
								vChoiceOne[j], vChoiceOne[jprime],
	                            vChoiceTwo[l], vChoiceTwo[lprime])
			# @error "to do"
		elseif rewardmat == :nobuild
			# need to be VERY careful with order of state vars here.. could get fucked up..
			reward = rewardfunc(p, getindex.(p.tStateVectors, [j, l, ixI...]),
								[vChoiceOne[jprime], vChoiceTwo[lprime]])
			# @error "to do"
		elseif rewardmat == :prebuild
			# jprime is first choice var, changes faster
			reward = mReward[jprime + nChoiceOne * (lprime-1), j + nChoiceOne * (l-1) + (nChoiceOne*nChoiceTwo) * (i-1)] # nChoices x nStates
		end

        if intdim == separable # mβEV is nChoices x nStochStates
            βEV = mβEV[jprime + nChoiceOne * (lprime-1), i] # mβEV is already discounted, jprime = j for inactive
        elseif intdim == intermediate # mβEV is nChoices x nStates
            βEV = mβEV[jprime + nChoiceOne * (lprime-1), j + nChoiceOne * (l-1) + (nChoiceOne*nChoiceTwo) * (i-1)]
        end

        valueProvisionalTwo = reward + βEV # mβEV is already discounted

        if valueProvisionalTwo >= valueHighSoFarTwo
           valueHighSoFarTwo = valueProvisionalTwo
           iChoice2 = lprime
        elseif concavity[2]
            break
        end

    end # lprime

    return valueHighSoFarTwo, iChoice2
end


function solve(p::TwoChoiceVar, method::Type{T},
		mTransition::Array{Float64,2}, mReward::Union{Array{Float64,2}, Nothing},
		disp::Bool, rewardmat::Symbol, monotonicity::Vector{Bool}, concavity::Vector{Bool}) where
			T <: Union{separable, intermediate}

	tChoiceVectors = p.tStateVectors[p.bEndogStateVars]
    (nChoiceOne, nChoiceTwo) = length.(tChoiceVectors)
	nChoices = nChoiceOne * nChoiceTwo
    (vChoiceOne, vChoiceTwo) = tChoiceVectors

    # nStochStates = size(mTransition,2)

    nStates = prod(length.(p.tStateVectors))
	tOtherStates = p.tStateVectors[.!p.bEndogStateVars]
    nOtherStates = prod(length.(tOtherStates))

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
	mTransition_βT = p.params.β * transpose(mTransition)

    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

        # @inbounds
		i = 0
		for ix in CartesianIndices(length.(tOtherStates)) # other states
			i = i + 1

            if monotonicity[2]
                vChoice2Start[:] .= 1
            end

            for l = 1:nChoiceTwo

                # We start from previous choice (monotonicity of policy function)
                if monotonicity[1]
                   iChoice1Start = 1
                end

                for j = 1:nChoiceOne # first state

                    valueHighSoFarOne = -Inf
                    valueProvisionalOne = -Inf
                    iChoice1 = 0
                    iChoice2 = 0

                    for jprime = iChoice1Start:nChoiceOne

                        (valueHighSoFarTwo, lprime) =
                                        getsecondchoice(p,
                                                        # active, # whether capital choice active or passive
                                                        mReward, mβEV, i, j, l,
                                                        nChoiceOne, vChoiceOne, jprime,
                                                        nChoiceTwo, vChoiceTwo,
                                                        vChoice2Start[j],
														ix.I,
														method,
														rewardmat,
														concavity)

                        valueProvisionalOne = valueHighSoFarTwo

                        if (valueProvisionalOne>=valueHighSoFarOne)
                            valueHighSoFarOne = valueProvisionalOne
                            iChoice1 = jprime
                            iChoice2 = lprime
                            if monotonicity[1]
                	          iChoice1Start = jprime
                            end
                        elseif concavity[1]
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

                    if monotonicity[2]
                        vChoice2Start[j] = iChoice2
                    end

                end #j (capital)

            end # l (second choice)

        end #i (stochastic variables)

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
            println("Parameters used:")
			println(p.params)
            break
        end

    end #while


    # println("Iteration ", string(iteration), ": Sup change is = ", string(maxDifference))

    mPolFun1 = vChoiceOne[mPolFunInd1]
    mPolFun2 = vChoiceTwo[mPolFunInd2]

    nNodes = length.(p.tStateVectors)
    meshPolFun1 = reshape(mPolFun1, tuple(nNodes...))
    meshPolFun2 = reshape(mPolFun2, tuple(nNodes...))

    meshValFun = reshape(mValFun, tuple(nNodes...))
    # meshExit   = reshape(mExit, tuple(nNodes...))

	return meshValFun, (meshPolFun1, meshPolFun2)
end # solve
