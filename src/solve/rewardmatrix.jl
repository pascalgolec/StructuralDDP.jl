function rewardmatrix(p::SingleChoiceVar)
    # output dimension is nChoices x nStates

    nChoices = length(p.vChoiceVector)

    nStates  = prod(length.(p.tStateVectors))
    mStates = gridmake(p.tStateVectors...)

    mF = zeros(Float64, nChoices, nStates)
    for i = 1 : nStates # state i
           for j = 1 : nChoices # choice j
               mF[j,i] = rewardfunc(p, mStates[i,:], p.vChoiceVector[j])
           end
    end

    return mF

end # getF
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
