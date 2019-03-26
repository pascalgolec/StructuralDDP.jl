rewardmatrix(p::DDM) =
    _rewardmatrix(p.rewardfunc, p.tStateVectors, p.tChoiceVectors, p.bEndogStateVars)

function _rewardmatrix(rewardfunc::Function, tStateVectors, tChoiceVectors, bEndogStateVars)
    # output dimension is nChoices x nStates

    tChoiceVectors = tStateVectors[bEndogStateVars]
    # nChoices = prod(length.(tChoiceVectors))
    mChoices = gridmake(tChoiceVectors...)
    nChoices = size(mChoices, 1)

    nStates  = prod(length.(tStateVectors))
    mStates = gridmake(tStateVectors...)

    mReward = zeros(Float64, nChoices, nStates)
    for i = 1 : nStates # state i
           for j = 1 : nChoices # choice j
               mReward[j,i] = rewardfunc(mStates[i,:], mChoices[j,:])
           end
    end

    return mReward
end

rewardmatrix_partial(p::DDM) = _rewardmatrix_partial(p.grossprofits,
    p.tStateVectors, p.bEndogStateVars)

function _rewardmatrix_partial(grossprofits::Function, tStateVectors, bEndogStateVars)
    # output dimension is nEndogStates x nExogStates

    # IMPORTANT THAT ENDOGENOUS STATE VARIABLES COME FIRST IN tStateVectors
    nEndogStates = prod(length.(tStateVectors[bEndogStateVars]))
    nExogStates = prod(length.(tStateVectors[.!bEndogStateVars]))

    mStates = gridmake(tStateVectors...)

    mReward = zeros(Float64, nEndogStates, nExogStates)
    for i = 1 : nExogStates # state i
           for j = 1 : nEndogStates # choice j
               # PROBABLY NEED TO CHANGE THIS, USE SOME CARTESIAN INDEXER
               mReward[j,i] = grossprofits(mStates[j + nEndogStates *(i-1),:])
           end
    end

    return mReward

end
