
# this function checks and dispatches
function solve(p::DDM; mTransition::Union{Nothing, Array{Float64,2}} = nothing,
    disp::Bool = false,
    rewardmat::Symbol = :nobuild,
    intdim::Symbol = :separable,
    monotonicity::Union{Bool,Vector{Bool}} = false,
    concavity::Union{Bool,Vector{Bool}} = false)

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
        # make both choices true or false if the user supplied only a single value or none
        if typeof(p)<:TwoChoiceVar
            if typeof(monotonicity) == Bool; monotonicity = fill!(Vector{Bool}(undef, 2), monotonicity) end
            if typeof(concavity) == Bool; concavity = fill!(Vector{Bool}(undef, 2), concavity) end
        end
        solve(p, eval(intdim), mTransition, mReward, disp, rewardmat, monotonicity, concavity)
    end
end
