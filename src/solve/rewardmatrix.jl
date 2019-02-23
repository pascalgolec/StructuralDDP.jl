function rewardmatrix(p::SingleChoiceVar)
    # output dimension is nChoices x nStates

    nChoices = length(p.tChoiceVectors[1])

    nStates  = prod(length.(p.tStateVectors))
    mStates = gridmake(p.tStateVectors...)

    mReward = zeros(Float64, nChoices, nStates)
    for i = 1 : nStates # state i
           for j = 1 : nChoices # choice j
               mReward[j,i] = rewardfunc(p, mStates[i,:], p.tChoiceVectors[1][j])
           end
    end

    return mReward

end

function rewardmatrix_partial(p::SingleChoiceVar)
    # output dimension is nEndogStates x nExogStates

    # IMPORTANT THAT ENDOGENOUS STATE VARIABLES COME FIRST IN tStateVectors

    nChoices = length(p.tChoiceVectors[1])

    nEndogStates = length(p.tStateVectors[ptest.bEndogStateVars][1])
    nExogStates = prod(length.(p.tStateVectors[.!ptest.bEndogStateVars]))

    mStates = gridmake(p.tStateVectors...)

    mReward = zeros(Float64, nEndogStates, nExogStates)
    for i = 1 : nExogStates # state i
           for j = 1 : nEndogStates # choice j
               mReward[j,i] = grossprofits(p, mStates[j + nEndogStates *(i-1),:])
           end
    end

    return mReward

end
#
# function rewardmatrix(p::ManyChoiceVars)
#
#
#     mChoices = gridmake(p.tChoiceVectors...)
#     (nChoices, nChoiceVars) = size(mChoices)
#
#     nStates  = prod(length.(p.tStateVectors))
#     mStates = gridmake(p.tStateVectors...)
#
#     mF = zeros(Float64, nChoices, nStates)
#     for i = 1 : nStates # state i
#            for j = 1 : nChoices # choice j
#                mF[j,i] = rewardfunc(p, mStates[i,:], mChoices[j,:])
#            end
#     end
#
#     return mF
#
# end # getF
