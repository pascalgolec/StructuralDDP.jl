
"""
Solve the dynamic programming problem.

##### Parameters
- `p::DDP` : Object that contains the dynamic optimization problem. See `createDDP` for details.
- `;disp::Bool(false)`: turn on/off displaying of intermeditate value function iteration steps
- `;disp_each_iter::Int(10)`: how often to display output
- `;max_iter::Int(500)`: Maximum number of iterations
- `;epsilon::Float64(1e-3)`: Value for epsilon-optimality. A lower value means the solution will be more accurate.
- `;rewardcall::Symbol(:jit)`:  Symbol specifying the handling of the reward in the VFI.
	Acceptable arguments are `:jit` if the reward matrix should not be prebuilt,
	`:pre` if it should be calculated beforehand, and `:pre_partial` if
	only part of the reward matrix that only depends on states should be calculated beforehand.
- `;monotonicity::Union{Bool,Vector{Bool}}(false)`: Option that can speed-up the solver
if optimal actions are monotonically increasing in their respective state variable.
Assumes that the action is next period's state.
- `;concavity::Union{Bool,Vector{Bool}}(false)`:  exploit the concavity of the value
function in the choice of next period's state.
- `:numquadnodes::Vector{Int}`: number of quadrature nodes to use for calculating expectations.

##### Returns

TODO

"""
function solve(p::DDP;
    disp::Bool = false,
    disp_each_iter::Int = 10,
    max_iter::Int = 500,
    epsilon::Float64 = 1e-3,
    rewardcall::Symbol = :jit,
    monotonicity::Union{Bool,Vector{Bool}} = false,
    concavity::Union{Bool,Vector{Bool}} = false,
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ)),
	mTransition::Union{Nothing, Array{Float64,2}} = nothing,
    intdim::Type{T} = p.intdim, # if want to override type in the model, mostly for testing
    ) where T<:IntDim

    if rewardcall == :pre || intdim == All
        mReward = rewardmatrix(p)
    elseif rewardcall == :pre_partial
        mReward = rewardmatrix_partial(p)
    elseif rewardcall == :jit
        mReward = nothing
    else
        error("supplied wrong rewardcall option")
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

        meshValFun, tmeshPolFun = _solve(p, intdim, mTransition, mReward,
            disp, disp_each_iter, max_iter, epsilon,
            rewardcall, monotonicity, concavity)
    end

    createsolution(p, meshValFun, tmeshPolFun)
end

function display_iter(iter::Int, disp_each_iter::Int, max_diff::Float64)
	if mod(iter, disp_each_iter)==0 || iter == 1
		println(" Iteration = ", iter, " Sup Diff = ", round(max_diff, sigdigits=5))
	end
end
