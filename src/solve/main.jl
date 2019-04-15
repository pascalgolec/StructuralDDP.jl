
"""
Solve the dynamic programming problem.

##### Parameters
- `p::DDM` : Object that contains the dynamic optimization problem. See `createDDP` for details.
- `;disp::Bool(false)`: turn on/off displaying of intermeditate value function iteration steps
- `;disp_each_iter::Int(10)`: how often to display output
- `;max_iter::Int(500)`: Maximum number of iterations
- `;epsilon::Float64(1e-3)`: Value for epsilon-optimality. A lower value means the solution will be more accurate.
- `;rewardmat::Symbol(:nobuild)`:  Symbol specifying the handling of the reward in the VFI.
	Acceptable arguments are `:nobuild` if the reward matrix should not be prebuilt,
	`:prebuild` if it should be calculated beforehand, and `:prebuild_partial` if
	only part of the reward matrix that only depends on states should be calculated beforehand.
- `;monotonicity::Union{Bool,Vector{Bool}}(false)`: Option that can speed-up the solver
if optimal actions are monotonically increasing in their respective state variable.
Assumes that the action is next period's state.
- `;concavity::Union{Bool,Vector{Bool}}(false)`:  exploit the concavity of the value
function in the choice of next period's state.
- `:initialize_exact::Bool(false)`: TODO
- `:numquadnodes::Vector{Int}`: number of quadrature nodes to use for calculating expectations.

##### Returns

TODO

"""
function solve(p::DDM;
    disp::Bool = false,
    disp_each_iter::Int = 10,
    max_iter::Int = 500,
    epsilon::Float64 = 1e-3,
    rewardmat::Symbol = :nobuild,
    monotonicity::Union{Bool,Vector{Bool}} = false,
    concavity::Union{Bool,Vector{Bool}} = false,
    initialize_exact::Bool = false,
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ)),
	mTransition::Union{Nothing, Array{Float64,2}} = nothing,
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

        meshValFun, tmeshPolFun = _solve(p, intdim, mTransition, mReward,
            disp, disp_each_iter, max_iter, epsilon,
            rewardmat, monotonicity, concavity)
    end

    createsolution(p, meshValFun, tmeshPolFun, initialize_exact)
end

function display_iter(iter::Int, disp_each_iter::Int, max_diff::Float64)
	if mod(iter, disp_each_iter)==0 || iter == 1
		println(" Iteration = ", iter, " Sup Diff = ", round(max_diff, sigdigits=5))
	end
end
