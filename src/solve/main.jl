
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

    if rewardmat == :prebuild || intdim == All
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

    if intdim == All
        meshValFun, tmeshPolFun = _solve(p, eval(intdim), mTransition, mReward, disp)
    else
        if typeof(p) <: DDP{NS,2} where NS
            # if user only supplied monotonivity/concavity for one option, extend to vector
            if typeof(monotonicity) == Bool; monotonicity = fill!(Vector{Bool}(undef, 2), monotonicity) end
            if typeof(concavity) == Bool; concavity = fill!(Vector{Bool}(undef, 2), concavity) end
        end

        meshValFun, tmeshPolFun = _solve(p, intdim, mTransition, mReward, disp,
            rewardmat, monotonicity, concavity)
    end

    createsolution(p, meshValFun, tmeshPolFun, initialize_exact)
end
