
# this function checks and dispatches
function solve(p::DDM; mTransition::Union{Nothing, Array{Float64,2}} = nothing,
    disp::Bool = false,
    rewardmat::Symbol = :nobuild,
    intdim::Symbol = :separable,
    monotonicity::Bool = false,
    concavity::Bool = false)

    intdim in(:SA, :intermediate, :separable) || error("supplied wrong intdim")

    if rewardmat == :prebuild || intdim == :SA
        mReward = rewardmatrix(p)
    elseif rewardmat == :prebuild_partial
        mReward = rewardmatrix_partial(p)
    elseif rewardmat == :nobuild
        mReward = nothing
    else
        error("supplied wrong rewardmat option")
    end

    if mTransition == nothing
        mTransition = Array(transitionmatrix(p, intdim = intdim))
    end

    if intdim == :SA
        solve(p, eval(intdim), mTransition, mReward, disp)
    else
        solve(p, eval(intdim), mTransition, mReward, disp, rewardmat, intdim, monotonicity, concavity)
    end
end


#
# function solve(p::DDM; mTransition::Union{Nothing, Array{Float64,2}} = nothing,
#     disp::Bool = false, intdim::Symbol=:separable)
#
#     @unpack rewardmat = p.params
#     if rewardmat == :prebuild_partial
#         mReward = outputfunc(p)
#     elseif rewardmat == :nobuild
#         mReward = nothing
#     elseif rewardmat == :prebuild
#         mReward = rewardmatrix(p)
#     end
#
#     @show intdim
#     if mTransition == nothing
#         mTransition = Array(transitionmatrix(p, intdim = intdim))
#     end
#     @show size(mTransition)
#
#     # not sure what the rule should be here...
#     # nOtherStates = prod(length.(p.tStateVectors[.!p.bEndogStateVars]))
#     # nFirstChoiceStates = length(p.tStateVectors[1])
#     #
#     # one = size(mOutput,1)
#     # two = nOtherStates
#     # one == two  ||
#     #     error("second dim of output matrix has wrong dimensions, is $one, should be $two")
#
#     out = solve(p, p.params.β * mTransition;
#         mReward = mReward, disp = disp, intdim = intdim)
# end
