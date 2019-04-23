# DiscreteDynamicProgramming.jl Documentation

This package solves and simulates discrete dynamic choice models with value function iteration. It has a simple and transparent syntax for defining dynamic optimization problems and is optimized for speed. The available options take advantage of properties of the problem which lead to orders of a magnitude speedup without relying on parallelization.

The general workflow is to define a problem, solve it, simulate it, and then analyze it. The full code for solving, simulating and plotting an example is:

```julia
using DiscreteDynamicProgramming
prob = CapitalAdjustModel(nK=150, nz=15, ρ=0.5, σ=0.3, γ=2.)
sol = solve(prob)
sim = simulate(sol, nPeriods=50, nFirms=5)
df_sim = DataFrame(sim)
using Plots, StatPlots
@df df_sim plot(:period, :state_1, group=:firm,
    xlabel="period", ylabel="capital stock")
```

where the individual pieces are described below.

## Step 1: define problem

Each infinite-horizon dynamic maximization problem has an objective function as follows

```math
\mathbb{E} \sum_{t=0}^{\infty} \beta^t r(s_t, a_t)
```

where

- ``s_t`` are the state variables and ``a_t`` the choice variable(s).
- ``r(s_t, a_t)`` is the current reward
- ``β`` is the discount factor

The state variables evolve according to a transition function ``\Phi(s, a, \varepsilon')`` depending on ``(s_t, a_t)`` and i.i.d shocks ``\varepsilon_{t+1}``. For a more formal definiton, see the [QuantEcon lectures](https://lectures.quantecon.org/jl/discrete_dp.html#Discrete-DPs). We can write the problem in a recursive form:

```math
V(s) = \max_a r(s, a) + \beta \mathbb{E} V(s') \\
\text{where } s' = \Phi(s, a, \varepsilon')
```

where ``s'`` corresponds to next period's state ``s_{t+1}``.

We will solve a capital investment model with convex adjustment costs as an example. It is defined as follows:

```math
V(K,z) = \max_I K^\alpha e^z - I - \frac{\gamma}{2} \frac{I^2}{K} + \beta \mathbb{E} V(K', z') \\
\text{where } K' = (1-\delta)K + I \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

where ``K`` is current and ``K'`` next period's capital, ``I`` is capital investment. Capital depreciates at the rate ``\delta``. Investment ``I`` adds to the capital stock but takes a period to become productive. The parameter ``\gamma`` determines the adjustment costs for investment, over and above the cost of capital ``I``. The productivity of capital is determined by ``z``, which follows an AR(1) process with autocorrelation ``\rho`` and volatility ``\sigma``. There are decreasing returns to scale for ``\alpha < 1``, so there is an optimal scale of the firm depending on ``z``. So the problem of the firm is to choose the optimal scale for next period, taking into account that next period's productivity ``z'`` is not perfectly known and large adjustments are disproportianally costly.

We define the problem by providing the state and action space, the reward function, the transition function and the factor at which future rewards are discounted.

```julia
β = 0.9 # discount factor
α = 0.67 # returns to scale parameter
ρ = 0.6 # persistence of producitivity
σ = 0.3 # volatility of productivity
δ = 0.15 # deprectiation rate of capital
γ = 2.0 # convex adjustment cost parameter


# define state and action space
nz = 5 # number of nodes for z
stdz = sqrt(σ^2/(1-ρ^2))
minz = -3*stdz
maxz =  3*stdz
vz = collect(LinRange(minz, maxz, nz))

nK = 100 # number of nodes for K
K_ss = α * exp(0.)/ ((1-β)/β*(1+γ*δ) + δ)) ^ (1/(1-α) # analytical steady state of K
log_K_min = 0.1 * log(K_ss)
log_K_max = 2 * log(K_ss)
vK   = exp.(collect(LinRange(log_K_min, log_K_max, nK))) # log-spaced grid for K

nI = 100 # number of nodes for I
K_diff = K_max - K_min
vI = collect(LinRange(-K_diff/2, K_diff/2, nI))

tStateVectors = (vK, vz)
tChoiceVectors = (vI,)


# define reward and transition functions
function reward(vStates, vChoices)
    K, z = vStates # note how the order corresponds to tStateVectors
    Kprime = vChoices[1]
    capx = Kprime - (1-δ)*K
    return K^α * exp(z) + γ/2*capx^2/K
end

function transition(vStates, vChoices, vShocks)
    K, z = vStates # note how the order corresponds to tStateVectors
    Kprime = vChoices[1]
    zprime  = ρ*z + σ * vShocks[1]
    return K, zprime
end

shockdistribution = Normal() # one-dimensional standard normal shock


# construct problem
prob = DiscreteDynamicProblem(
            tStateVectors,
            tChoiceVectors,
            reward,
            transition,
            shockdistribution,
            β)
```

The solver supports any type of univariate distribution distribution and multivariate-normal distributions. The [Distributions.jl package](https://juliastats.github.io/Distributions.jl/stable/) contains the syntax for implementing these distributions.

The state space is strictly enforced - the state variables are always forced to stay within their bounds upper and lower bounds.

Note that in the problem definition above, the transition of the state variables is a function of states, choices and shocks. The solver will integrate over a high dimension of variables when calculating expectations, which is rather slow. The [problem options section](#Problem-Options-1) explains how to reformulate a large family of problems including the one above to gain orders of magnitude speedup.

## Step 2: solve problem

We can solve the model with the function `solve`, which returns the (discrete) value and optimal policy for each point in the state-space.

```julia
sol = solve(prob)
```

The result of `solve` is a solution object. We can access the value and policy of the fifth grid point of the first state and the third grid point of the second state by:
```julia
value(sol)[5, 3]
policy(sol)[5, 3]
```

If there are multiple choice variables then `policy(sol)` returns a tuple, so we would have to specify that we for example are interested in the second choice variable:
```julia
policy(sol)[2][5, 3]
```

The `policy` and `value` objects that are returned by default act as a continuous solution via an interpolation. We can access the interpolated values by treating them as functions, for example:
```julia
policy(sol)(10., 0.5) # the optimal policy for state one = 10. and state two = 0.5
```
Note the difference between these: indexing with `[i]` is the policy at the ith grid point, while `(k)` is an interpolation for state `k`. Also note that if an interpolation outside of the state grid is requested, then the value/policy at the closest grid points is returned instead.

The solver can be controlled using different options which are discribed in the [Solver Options section](#Solver-Options-1). For example, we can tell the solver to precompute the reward for the different combinations of states and choices before starting the value function iteration:

```julia
sol = solve(prob; rewardmat=:prebuild)
```

Experimenting with the solver options is essential if one wants to solve the model as fast as possible.

## Step 3: simulate

We can simulate a panel:

```julia
sim = simulate(sol, nPeriods = 60, nFirms = 100)
```

The function `simulate` returns a simulation object, from which we can retreive the simulated states, choices and value as individual arrays:

```julia
sim_states = states(sim)
sim_policy = policy(sim)
sim_val = value(sim)
```

Convenience methods which convert the simulation object into an array or dataframe:
```julia
a_sim = Array(sim)
df_sim = DataFrame(sim)
```

We can then use the dataframe to analyse the simulation, for example by plotting a histogram:
```julia
using Plots
histogram(df_sim.state_1, xlabel="K",legend=false)
```

Or the time of capital of each firm takes:
```julia
using StatPlots
@df df_sim plot(:period, :state_1, group=:firm, xlabel="time", ylabel="K")
```

By default, all firms start with the same state variables in the simulation. More details about this and more sophisticated starting points are in the [simulator options section](#Simulator-Options-1).

# Problem Options

## Lower integration dimension for speed

The solver must calculate expectations of future states. By default, it must integrate over all state variables, choice variables and shocks, i.e. ``s_{t+1} = \Phi(s_t, a_t, \varepsilon_{t+1})``. Many problems can be rewritten in a way such that we can integrate over fewer variables, which leads to orders of magnitude speedup.

The first step is to rewrite the problem such that the action is equal to the corresponding next period's state. In our example, the action ``a`` of the firm then is not investment but simply next period's capital stock ``K'``.

```math
V(K,z) = \max_a K^\alpha e^z - (a - (1-\delta) K) - \frac{\gamma}{2} \frac{(a - (1-\delta) K)^2}{K} + \beta \mathbb{E} V(K', z') \\
\text{where } K' = a \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

By separating endogenous and exogenous state variables, resp. ``K`` and ``z``, this way, we can define the problem as follows:

```julia
tStateVectors = (vK, vz)
tChoiceVectors = (1,) # specify which state corresponds to the choice
function transfunc(vStates, vChoices, vShocks)
    z = vStateVars[2]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime # only integrate z
end
prob = DiscreteDynamicProblem(...,
            transfunc,
            tChoiceVectors;
            intdim = :Separable) # specify the integration dimension
```

We need to specify which state corresponds to the choice, only return ``z'`` in the transition function and specify that our integration dimension is of the type `:Separable`.

In our example however we can do even better if we tell the solver that the choice of `Kprime` does not affect the transition of ``z``.

```julia
tStateVectors = (vK, vz)
tChoiceVectors = (1,)
function transfunc(vStates, vShocks) # no choices as input
    z = vStateVars[2]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
prob = DiscreteDynamicProblem(...,
            transfunc,
            tChoiceVectors;
            intdim = :Separable_States) # specify the integration dimension
```

Note that the transition function only takes states and shocks as an input and that the integration dimension now is of the type `:Separable_States`. In our example however we can do even better if we tell the solver that also the endogenous state ``K`` does not affect the transition of ``z``.

```julia
function transfunc(vExogStates, vShocks)
    z = vExogStates[1] # note the different index for z
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
prob = DiscreteDynamicProblem(...,
            transfunc,
            tChoiceVectors;
            intdim = :Separable_ExogStates) # specify the integration dimension
```

Note that the transition function only takes the exogenous states and shocks as an input and that the integration dimension now is of the type `:Separable_ExogStates`.


The following summarises the four different types of integration dimensions:

- `:All` - next period states as a function of states, actions and shocks, i.e. ``s' = \Phi(s, a, \varepsilon')``. This formulation is the slowest and corresponds to the baseline example.
- `:Separable` - next period exogenous states as a function of all states, actions and shocks, i.e. ``s'_e = \Phi(s, a, \varepsilon')``.
- `:Separable_States` - next period exogenous states as a function of all states and shocks, i.e. ``s'_e = \Phi(s, \varepsilon')``.
- `:Separable_ExogStates` - next period exogenous states as a function of only exogenous states and shocks, i.e. ``s'_e = \Phi(s_e, \varepsilon')``.

For the solver to be fast one should use the lowest possible dimension. For all separable integration dimensions one must specify which state vector corresponds to the choice vector, i.e. `tChoiceVectors = (1,)`. Important: when defining the state space `tStateVectors`, the endogenous state variable(s) must come first.

## Concise notation

In some problems one may only have one shock or one choice variable. In this case, one can write the reward- and transition functions more concisely in terms as numbers as inputs as opposed to vectors:

```julia
function reward(vStateVars, Kprime)
    (K, z) = vStateVars
    capx = Kprime - (1-δ)*K
    return K^α * exp(z) + γ/2*(capx/K- δ)^2 * K
end
function transition(z, shock)
    return zprime  = ρ*z + σ * shock
end
```

# Solver Options

## For speed

There are different options available which make solving the model faster. One should check whether one can lower the integration dimension of the problem first though, described in [the problem options section](#Problem-Options-1).

### Monotonicity and concavity

The `monotonicity` keyword allows the solver to exploit the monotonicity of the choice of next period's state in the current state. Put differently, if ``s^{d*}_{t+1}`` is the optimal choice for ``s^{d*}_t``, then the optimal choice for ``s^{d'}_t \geq s^{d*}_t`` is ``s^{d'}_{t+1} \geq s^{d*}_{t+1}`` (keeping fixed the other state variables). The default is `false`.

The `concavity` keyword allows the solver to exploit the concavity of the value function in the choice of next period's state. Put differently, if ``V(s^{d'}_{t+1}, s^d_t, s^s_t) < V(s^{d*}_{t+1}, s^d_t, s^s_t)``, then also ``V(s^d_{t+1}, s^d_t, s^s_t) < V(s^{d'}_{t+1}, s^d_t, s^s_t)`` for any ``s^d_{t+1} > s^{d'}_{t+1}``. The default is `false`.

A handy way to check whether the problem fulfils these conditions is the `compare` function. It checks whether tow different solutions are identical:

```julia
sol = solve(prob; concavity=false)
sol_conc = solve(prob; concavity=true)
compare(sol, sol_conc; tol=1e-4)
```

Note: `monotonicity` and `concavity` currently only work if the integration dimension is separable. The options don't work if there are discontinuities in the reward function.


### Prebuilding the reward matrix

The `rewardcall` keyword determines whether (part of) the reward for each state-action pair should be precomputed before the value function iteration. Prebuilding requires more memory but less computations during the VFI.

`:jit` - the reward is calculated during the VFI whenever it is required, i.e. just in time. This requires little memory but more computations during the VFI. This is the default.

`:pre` - precompute the reward for each possible combination of states and actions before the VFI. This requires more memory but less computations during the VFI.

`:pre_partial` - precompute part of the reward before the VFI that only depends on states, but not choices. From experience, this option is the fastest when combined with monotonicity and concavity. To exploit this option, we must supply the inner and outer part of the reward function when defining the problem. The outer part is the same argument in `DDP` as the standard reward function and the partial reward function enters as a keyword arguement `rewardfunc_partial`. In the solver we must then specify that the reward should be partially precomputed.

```julia
function reward(Partial_Reward, vStateVars, Kprime)
    K = vStateVars[1]
    capx = Kprime - (1-δ)K
    return (1-τ)*(Partial_Reward - F*K - γ/2*(capx/K- δ)^2 * K) - (1-κ*(capx<0))*capx + τ * δ * K
end
reward_partial(vStateVars) = vStateVars[1]^α * exp(vStateVars[2])

prob = DiscreteDynamicProblem(
    tStateVectors,
    tChoiceVectors,
    reward,
    transition,
    shockdistribution,
    β;
    rewardfunc_partial = reward_partial)

solve(prob, rewardcall = :pre_partial)
```


## Transition matrix

The solver calculates a markov transition matrix from the transition function using quadrature methods. We can control the number of quadrature nodes for each shock as with the keyword argument `numquadnodes`. The default is `[5]`.

It is also possible to supply your own transition matrix with the keyword argument `mTransition`:

```julia
using QuantEcon: tauchen
sol = solve(prob, mTransition = tauchen(nz, ρ, σ))
```

## Miscellaneous

- `epsilon`: Value for epsilon-optimality. Determines how accurate the solution is. Default is 1e-3.
- `max_iter`: Maximum number of iterations before stopping. Defaults to 1e5.
- `disp`: whether to display number of iterations and epsilon-convergence during solving. Default is false.
- `disp_each_iter`: wait how many iterations for displaying status during solving. Default is 10.

# Simulator Options

The keyword `get_value` determines whether the simulator should also compute the value at the beginning of each time period. This takes more time. The default is `get_value = false.`.

## Initialization

The default setting how the initial state variables are drawn at t=0 is the burn-in method. Here, each firm starts with states that are the median values of the state space. The idea is to simulate the model for more periods than necessary and then only analyse the simulated values after 60 periods or so on, depending on the model.

There is also an option for more sophisticated initialization, where the initial endogenous state variables are chosen by the firm subject to a loss function and the exogenous ones are predetermined or random. For the neoclassical model, this entails coding the following functions as an additional input into the problem:

```julia
initprob(value, vChoices) = value - (1 + (1-β)/β + C0) * vChoices[1]
init(vShocks) = vShocks[1] * sqrt(σ^2 / (1-ρ^2)) # = z0
createDiscreteDynamicProblem(<other variables>;
    initializationproblem = initprob,
    initializefunc = init,
    tChoiceVectorsZero = (1,))
```

If we have a problem where our choice variables are not exacly equal to some of the state variables, then we need to specify additionally which state variables the firm chooses at t=0. This is achieved with the option `tChoiceVectorsZero`, i.e. `tChoiceVectorsZero = (1,)` for the neoclassical model. If `tChoiceVectorsZero` is not specified, then the solver/simulator uses indices of the choice variables for the dynamic optimization.

## Fixed shock draws

It is also possible to simulate a model solution using the identical draw of shocks:

```julia
shocks = drawshocks(prob, nPeriods = 60, nFirms = 100)
sim = simulate(sol, shocks)
```

Note: the draw of shocks refers to the supplied shock distribution in the problem defintion. If the shock distribution is parametrized, for example by its variance, then one should not do comparative statics on those parameters.

# notes

## To do

- [ ] support more than two choice variables when the integration dimension is `separable`
- [ ] plotrecipe for `DDPSolution`
- [ ] allow parametric types in problem definition of state and choice vectors, i.e. `AbstractVector{T}` instead of `Vector{Float64}`
- [ ] allow `LabelledArrays.jl` in problem definition
- [ ] allow discontinuities in the reward function with `monotonicity` and `concavity`
