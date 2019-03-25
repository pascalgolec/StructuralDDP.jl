function rewardmatrix(p::DDM)
    # output dimension is nChoices x nStates

    tChoiceVectors = p.tStateVectors[p.bEndogStateVars]
    # nChoices = prod(length.(tChoiceVectors))
    mChoices = gridmake(tChoiceVectors...)
    nChoices = size(mChoices, 1)

    nStates  = prod(length.(p.tStateVectors))
    mStates = gridmake(p.tStateVectors...)

    mReward = zeros(Float64, nChoices, nStates)
    for i = 1 : nStates # state i
           for j = 1 : nChoices # choice j
               mReward[j,i] = p.rewardfunc(mStates[i,:], mChoices[j,:])
           end
    end

    return mReward
end


function rewardmatrix_partial(p::DDM)
    # output dimension is nEndogStates x nExogStates

    # IMPORTANT THAT ENDOGENOUS STATE VARIABLES COME FIRST IN tStateVectors
    nEndogStates = prod(length.(p.tStateVectors[p.bEndogStateVars]))
    nExogStates = prod(length.(p.tStateVectors[.!p.bEndogStateVars]))

    mStates = gridmake(p.tStateVectors...)

    mReward = zeros(Float64, nEndogStates, nExogStates)
    for i = 1 : nExogStates # state i
           for j = 1 : nEndogStates # choice j
               # PROBABLY NEED TO CHANGE THIS, USE SOME CARTESIAN INDEXER
               mReward[j,i] = p.grossprofits(mStates[j + nEndogStates *(i-1),:])
           end
    end

    return mReward

end
