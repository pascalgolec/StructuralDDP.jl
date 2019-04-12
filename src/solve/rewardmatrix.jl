rewardmatrix(p::DDP{nStateVars,nChoiceVars,C}) where
    {nStateVars,nChoiceVars,C<:Int} =
    _rewardmatrix(p.rewardfunc, p.tStateVectors, getchoicevars(p.tStateVectors, p.tChoiceVectors))

rewardmatrix(p::DDP{nStateVars,nChoiceVars,C}) where
    {nStateVars,nChoiceVars,C<:Vector{Float64}} =
    _rewardmatrix(p.rewardfunc, p.tStateVectors, p.tChoiceVectors)

""" Return a nChoices x nStates reward matrix"""
function _rewardmatrix(rewardfunc::Function, tStateVectors, tChoiceVectors)
    # output dimension is nChoices x nStates

    nChoices = prod(length.(tChoiceVectors))
    nStates = prod(length.(tStateVectors))
    mReward = zeros(Float64, nChoices, nStates)

    iterator_states = Iterators.product(tStateVectors...)

    if length(tChoiceVectors) > 1
        iterator_choices = Iterators.product(tChoiceVectors...)
    else # don't want a tuple of choices if only one
        iterator_choices = tChoiceVectors[1]
    end

    for (i, states) in enumerate(iterator_states)
        for (j, choices) in enumerate(iterator_choices)
           # @show choices
           mReward[j,i] = rewardfunc(states, choices)
           # error("stop")
        end
    end

    return mReward
end

rewardmatrix_partial(p::DDM) = _rewardmatrix_partial(p.grossprofits,
    getchoicevars(p.tStateVectors, p.tChoiceVectors),
    getnonchoicevars(p.tStateVectors, p.tChoiceVectors))

"""Return a nEndogStates x nExogStates reward matrix."""
function _rewardmatrix_partial(grossprofits::Function,
    tEndogStateVectors, tExogStateVectors)
    # output dimension is nEndogStates x nExogStates

    # julia does not know dim of endog state vars from this expression
    nEndogStates::Int64 = prod(length.(tEndogStateVectors))
    nExogStates::Int64 = prod(length.(tExogStateVectors))

    # mStates::Array{Float64,2} = gridmake(tStateVectors...)

    mReward = zeros(Float64, nEndogStates, nExogStates)

    iterator_endogstates = Iterators.product(tEndogStateVectors...)
    iterator_exogstates = Iterators.product(tExogStateVectors...)

    for (i, exogstates) in enumerate(iterator_exogstates)
        for (j, endogstates) in enumerate(iterator_endogstates)
           mReward[j,i] = grossprofits((endogstates..., exogstates...))
        end
    end

    return mReward

end
