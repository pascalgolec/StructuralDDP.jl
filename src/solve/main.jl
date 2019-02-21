
function solve(p::DDM; mTransition::Union{Nothing, Array{Float64,2}} = nothing,
    disp::Bool = false)

    @unpack rewardmat = p.params
    if rewardmat == :prebuild_partial
        mReward = outputfunc(p)
    elseif rewardmat == :nobuild
        mReward = nothing
    elseif rewardmat == :prebuild
        mReward = rewardmatrix(p)
    end


    if mTransition == nothing
        mTransition = Array(transitionmatrix(p))
    end

    # not sure what the rule should be here...
    # nOtherStates = prod(length.(p.tStateVectors[.!p.bEndogStateVars]))
    # nFirstChoiceStates = length(p.tStateVectors[1])
    #
    # one = size(mOutput,1)
    # two = nOtherStates
    # one == two  ||
    #     error("second dim of output matrix has wrong dimensions, is $one, should be $two")

    out = solve(p, p.params.β * mTransition, mReward = mReward, disp = disp)
end
