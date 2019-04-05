rewardmatrix(p::DDP{nStateVars,nChoiceVars,C,G,IP,IF}) where
    {nStateVars,nChoiceVars,C<:Int64,G,IP,IF} =
    _rewardmatrix(p.rewardfunc, p.tStateVectors, getindex(p.tStateVectors, collect(p.tChoiceVectors)))
rewardmatrix(p::DDP{nStateVars,nChoiceVars,C,G,IP,IF}) where
    {nStateVars,nChoiceVars,C<:Vector{Real},G,IP,IF} =
    _rewardmatrix(p.rewardfunc, p.tStateVectors, p.tChoiceVectors)

function _rewardmatrix(rewardfunc::Function, tStateVectors, tChoiceVectors)
    # output dimension is nChoices x nStates

    mChoices = gridmake(tChoiceVectors...)::Array{Float64,2} # gridmake is not type stable
    nChoices = size(mChoices, 1)

    nStates  = prod(length.(tStateVectors))
    mStates = gridmake(tStateVectors...)

    mReward = zeros(Float64, nChoices, nStates)
    for i = 1 : nStates # state i
           for j = 1 : nChoices # chice j
               mReward[j,i] = rewardfunc(mStates[i,:], mChoices[j,:])
           end
    end

    return mReward
end

rewardmatrix_partial(p::DDM) = _rewardmatrix_partial(p.grossprofits,
    p.tStateVectors, p.bEndogStateVars)

function _rewardmatrix_partial(grossprofits::Function, tStateVectors, bEndogStateVars)
    # output dimension is nEndogStates x nExogStates

    # julia does not know dim of endog state vars from this expression
    nEndogStates::Int64 = prod(length.(tStateVectors[bEndogStateVars]))
    nExogStates::Int64 = prod(length.(tStateVectors[.!bEndogStateVars]))

    mStates::Array{Float64,2} = gridmake(tStateVectors...)

    mReward = zeros(Float64, nEndogStates, nExogStates)

    for i = 1 : nExogStates # state i
           for j = 1 : nEndogStates # choice j
               # COULD CHANGE THIS, USE SOME CARTESIAN INDEXER
               mReward[j,i] = grossprofits(mStates[j + nEndogStates *(i-1),:])
           end
    end

    return mReward

end
