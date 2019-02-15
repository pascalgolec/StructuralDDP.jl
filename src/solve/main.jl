
# Structure:
# - `p.bConvex==false` takes reward matrix as input, the others the output of the firm

solve(p::DDM; disp::Bool = false) = solve(p, eval(p.params.intdim), disp = disp)

function solve(p::DDM, method::Type{T}; disp::Bool = false) where
                        T <: Union{separable, intermediate}

    @unpack rewardmat = p.params
    if rewardmat == :prebuild_partial
        mReward = outputfunc(p)
    elseif rewardmat == :nobuild
        mReward = nothing
    elseif rewardmat == :prebuild
        mReward = rewardmatrix(p)
    end
    mG = transitionmatrix(p)

    # not sure what the rule should be here...
    # nOtherStates = prod(length.(p.tStateVectors[.!p.bEndogStateVars]))
    # nFirstChoiceStates = length(p.tStateVectors[1])
    #
    # one = size(mOutput,1)
    # two = nOtherStates
    # one == two  ||
    #     error("second dim of output matrix has wrong dimensions, is $one, should be $two")

    out = solve(p, p.params.β * Array(mG), mReward = mReward, disp = disp)
end
