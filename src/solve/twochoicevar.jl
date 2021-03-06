getreward(rewardcall::Type{pre_partial}, rewardfunc, mReward, tStateVectors,
	vChoiceOne, vChoiceTwo,
	cnt, ix, j, jprime, l, lprime) =
	rewardfunc(mReward[cnt.i_state],
					   getindex.(tStateVectors, (j, l, ix.I...)),
					   (vChoiceOne[jprime], vChoiceTwo[lprime]))

getreward(rewardcall::Type{jit}, rewardfunc, mReward, tStateVectors,
	vChoiceOne, vChoiceTwo,
	cnt, ix, j, jprime, l, lprime) =
	reward = rewardfunc(getindex.(tStateVectors, (j, l, ix.I...)),
						(vChoiceOne[jprime], vChoiceTwo[lprime]))

getreward(rewardcall::Type{pre}, rewardfunc, mReward, tStateVectors,
	vChoiceOne, vChoiceTwo,
	cnt, ix, j, jprime, l, lprime) =  mReward[cnt.i_choice, cnt.i_state]

_solve(p::DDP,
	mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
	opts::SolverOptions) =
	_solve(p, mTransition, mReward, opts, p.tStateVectors,
		getchoicevars(p.tStateVectors, p.tChoiceVectors),
		getnonchoicevars(p.tStateVectors, p.tChoiceVectors),)


function _solve(p::DDP,
				mTransition::Array{Float64,2}, mReward::Union{Array{Float64}, Nothing},
				opts::SolverOptions,
				tStateVectors,
				tChoiceVectors,
				tOtherStateVectors)

	@unpack β, rewardfunc = p
	@unpack disp, disp_each_iter, max_iter, epsilon, rewardcall,
	monotonicity, concavity, intdim = opts
	monotonicity1, monotonicity2 = monotonicity
	concavity1, concavity2 = concavity

    (nChoiceOne, nChoiceTwo) = length.(tChoiceVectors)
    (vChoiceOne, vChoiceTwo) = tChoiceVectors

	nChoices = nChoiceOne * nChoiceTwo
    nStates = prod(length.(tStateVectors))
    nOtherStates = prod(length.(tOtherStateVectors))

	mValFun    = zeros((nChoices, nOtherStates)) # matrix so that can multiply with mTransition
	mValFunNew = zeros(nStates) # vector so that can index easily

    mPolFunInd1 = zeros(Int16, nStates)
    mPolFunInd2 = zeros(Int16, nStates)

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

    vChoice2Start = ones(Int64, nChoiceOne)
    vChoice1Start = 1
	iChoice1inner = 0

	i::Int64 = 0 # outer loop for other state vars
	cnt = initialize_counter()

	# will need beta times transpose of transition matrix
	mTransition_βT = β * transpose(mTransition)

    # VFI
    while maxDifference > tolerance

        mul!(mβEV, mValFun, mTransition_βT)

		i = 0
		reset!(cnt)

		# @inbounds
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

					cnt.i_statechoice += (vChoice2Start[j]-1)*nChoiceOne # compensate if start later
					cnt.i_choice += (vChoice2Start[j]-1)*nChoiceOne # compensate if start later

                    for lprime = vChoice2Start[j]:nChoiceTwo

						# get optimal second choice variable conditional on first
						valueHighSoFarOne = -Inf
						valueProvisionalOne = -Inf
						iChoice1inner  = 0

						# find highest value for second state var
						cnt.i_statechoice += vChoice1Start-1 # compensate if start later
						cnt.i_choice += vChoice1Start-1
						for jprime = vChoice1Start:nChoiceOne

							cnt.i_statechoice += 1
							cnt.i_choice += 1

							reward = getreward(rewardcall, rewardfunc,
								mReward, tStateVectors,
								vChoiceOne, vChoiceTwo,
								cnt, ix, j, jprime, l, lprime)

							valueProvisionalOne = reward + mβEV[cnt.i_choice,
								getcounter(cnt, intdim)]

							if valueProvisionalOne >= valueHighSoFarOne
    							   valueHighSoFarOne = valueProvisionalOne
    							   iChoice1inner = jprime
	                        elseif concavity1
								adj = nChoiceOne - jprime # adjust counter if finish early
								cnt.i_statechoice += adj
								cnt.i_choice += adj
    							break
							end

						end # jprime

                        valueProvisionalTwo = valueHighSoFarOne

                        if valueProvisionalTwo > valueHighSoFarTwo
                            valueHighSoFarTwo = valueProvisionalTwo
                            iChoice1 = iChoice1inner
                            iChoice2 = lprime
                        elseif concavity2
							adj = (nChoiceTwo - lprime)*nChoiceOne # adjust counter if finish early
							cnt.i_statechoice += adj
							cnt.i_choice += adj
                            break
                        end

                    end #lprime

                    mValFunNew[cnt.i_state] = valueHighSoFarTwo
                    mPolFunInd1[cnt.i_state] = iChoice1
                    mPolFunInd2[cnt.i_state] = iChoice2

					if monotonicity1
                        vChoice1Start = iChoice1inner
                    end
					if monotonicity2
					  vChoice2Start[j] = iChoice2
					end

                end # j (first choice)

            end # l (second choice)

        end #i (stochastic variables)

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
