
# this function checks and dispatches
function solve(p::DDM; mTransition::Union{Nothing, Array{Float64,2}} = nothing,
    disp::Bool = false,
    rewardmat::Symbol = :nobuild,
    monotonicity::Union{Bool,Vector{Bool}} = false,
    concavity::Union{Bool,Vector{Bool}} = false,
    initialize_exact::Bool = false,
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ)),
    intdim::Type{T} = p.intdim, # if want to override type in the model, mostly for testing
    ) where T<:DDMIntDim

    # intdim in(:SA, :intermediate, :separable) || error("supplied wrong intdim")

    if rewardmat == :prebuild || intdim == SA
        mReward = rewardmatrix(p)
    elseif rewardmat == :prebuild_partial
        mReward = rewardmatrix_partial(p)
    elseif rewardmat == :nobuild
        mReward = nothing
    else
        error("supplied wrong rewardmat option")
    end

    if mTransition == nothing
        mTransition = Array(transitionmatrix(p, numquadnodes = numquadnodes))
    end

    if intdim == SA
        meshValFun, tmeshPolFun = _solve(p, eval(intdim), mTransition, mReward, disp)
    else
        # if typeof(p.tChoiceVectors) == Tuple{Vector{Float64}}
        #     meshValFun, tmeshPolFun = _solve1(p, eval(intdim), mTransition, mReward, disp,
        #         rewardmat, monotonicity, concavity)
        # elseif typeof(p.tChoiceVectors) == Tuple{Vector{Float64},Vector{Float64}}
        #
        #     # if user only supplied monotonivity/concavity for one option, extend to vector
        #     if typeof(monotonicity) == Bool; monotonicity = fill!(Vector{Bool}(undef, 2), monotonicity) end
        #     if typeof(concavity) == Bool; concavity = fill!(Vector{Bool}(undef, 2), concavity) end
        #
        #     meshValFun, tmeshPolFun = _solve2(p, eval(intdim), mTransition, mReward, disp,
        #         rewardmat, monotonicity, concavity)
        # else
        #     error("only :SA supports more than two choice variables")
        # end

        meshValFun, tmeshPolFun = _solve(p, eval(intdim), mTransition, mReward, disp,
            rewardmat, monotonicity, concavity)
    end

    createsolution(p, meshValFun, tmeshPolFun, initialize_exact)
end
