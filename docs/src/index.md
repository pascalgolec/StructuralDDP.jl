# StructuralDDP.jl Documentation

This package defines, solves and simulates discrete dynamic optimization problems with value function iteration. It features different options to add structure to the problem and provide information about its properties. A wide range of problems can be solved very fast this way.

The general workflow is to define a problem, solve it, simulate it, and then analyze it. The full code for an example is:

```julia
using StructuralDDP
prob = CapitalAdjustModel(nK=150, nz=15, ρ=0.5, σ=0.3, γ=2.)
sol = solve(prob)
sim = simulate(sol, nPeriods=50, nFirms=5)
df_sim = DataFrame(sim)
using Plots, StatPlots
@df df_sim plot(:period, :state_1, group=:firm,
    xlabel="years", ylabel="capital stock")
```

where the individual pieces are described below.

## Step 1: define problem

Each infinite-horizon dynamic maximization problem has an objective function as follows

```math
\mathbb{E} \sum_{t=0}^{\infty} \beta^t r(s_t, a_t)
```

where ``s_t`` are the state variables and ``a_t`` the choice variable(s), ``r(s_t, a_t)`` is the current reward and ``β`` is the discount factor.

The state variables evolve according to a transition function ``\Phi(s_t, a_t, \varepsilon_{t+1})`` depending on ``(s_t, a_t)`` and i.i.d shocks ``\varepsilon_{t+1}``. For a more formal definiton, see the [QuantEcon lectures](https://lectures.quantecon.org/jl/discrete_dp.html#Discrete-DPs). We can write the problem in a recursive form:

```math
V(s) = \max_a r(s, a) + \beta \mathbb{E} V(s') \\
\text{where } s' = \Phi(s, a, \varepsilon')
```

where ``s`` and ``s'`` correspond to the current state ``s_t`` and next period's state ``s_{t+1}``, respectively.

We will solve a capital investment model with convex adjustment costs as an example. It is defined as follows:

```math
V(K,z) = \max_i K^\alpha e^z - i K - \frac{\gamma}{2} i^2 K + \beta \mathbb{E} V(K', z') \\
\text{where } K' = (1-\delta + i) K \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

where ``K`` is current and ``K'`` next period's capital and ``i`` is the (net) capital investment rate. Capital depreciates at the rate ``\delta``. Investment ``i`` adds to the capital stock but takes a period to become productive. The parameter ``\gamma`` determines the adjustment costs for investment, over and above the cost of capital. The productivity of capital is determined by ``z``, which follows an AR(1) process with autocorrelation ``\rho`` and volatility ``\sigma``. There are decreasing returns to scale for ``\alpha < 1``, so there is an optimal scale of the firm depending on ``z``. The problem of the firm is to choose the optimal scale for next period.

We set up the problem by defining the state and action space, the reward function, the transition function and the factor at which future rewards are discounted.

```julia
using StructuralDDP, Distributions

# CapitalAdjustModel
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
K_ss = α * exp(0.)/ ((1-β)/β*(1+γ*δ) + δ) ^ (1/(1-α)) # analytical steady state of K
log_K_min = 0.1 * log(K_ss)
log_K_max = 2 * log(K_ss)
vK   = exp.(collect(LinRange(log_K_min, log_K_max, nK))) # log-spaced grid for K

ni = 100 # number of nodes for i
vi = collect(LinRange(-0.5, 2., ni))

tStateVectors = (vK, vz)
tChoiceVectors = (vi,)


# define reward and transition functions
function reward(vStates, vChoices)
    K, z = vStates # note how the order corresponds to tStateVectors
    i = vChoices[1]
    return K^α * exp(z) - i*K - γ/2 * i^2 * K
end

function transition(vStates, vChoices, vShocks)
    K, z = vStates # note how the order corresponds to tStateVectors
    i = vChoices[1]
    Kprime = (1-δ+i)*K
    zprime  = ρ*z + σ * vShocks[1]
    return Kprime, zprime
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

The `DiscreteDynamicProblem` constructor creates a problem instance. It is important that the order of the arguments in the `reward` function, `transition` function and `DiscreteDynamicProblem` are the same as in the example above.

The package supports any type of univariate shock distribution and multivariate-normal shock distributions. [Distributions.jl](https://juliastats.github.io/Distributions.jl/stable/) contains the syntax for implementing different shock distributions.

Note that in the problem definition above, the transition of the state variables is a function of states, choices and shocks. The solver will integrate over a high dimension of variables when calculating expectations, which is rather slow. The [problem options section](#Problem-Options-1) explains how to reformulate a large family of problems including the one above to gain orders of magnitude speedup.

The definition of the discrete state and action space means that the problem is solved for each point in the state space, where there finite number of possible actions. When simulating the model, linear interpolation is used to approximate the solution between grid points. Note that the state space is strict - the state variables are always forced to stay within their bounds upper and lower bounds.

## Step 2: solve problem

We can solve the model with the function `solve`, which returns the value and optimal policy for each point in the state-space.

```julia
sol = solve(prob)
```

The result of `solve` is a solution object. We can retreive the value function with `value` and the policy function with `policy`. For example, the value and policy for the fifth grid point of the `K` and the third grid point of `z` are:
```julia
value(sol)[5, 3] # 18.9
policy(sol)[5, 3] # 0.611
```

If there are multiple choice variables then `policy(sol)` returns a tuple, so we would have to specify that we for example are interested in the second choice variable:
```julia
policy(sol)[2][5, 3]
```

The `policy` and `value` objects that are returned by default act as a continuous solution via an interpolation. We can access the interpolated values by treating them as functions, for example:
```julia
policy(sol)(10.2, 0.5) # 0.273 = optimal policy for K = 10.2 and z = 0.5
```
Note the difference between these: indexing with `[i,j]` is the policy at the (i,j)th grid point, while `(K,z)` is an interpolation for state $(K,z)$. Also note that if an interpolation outside of the state grid is requested, then the value/policy at the closest grid points is returned instead.

The solver can be controlled using different options which are discribed in the [Solver Options section](#Solver-Options-1). For example, we can tell the solver to precompute the reward for all different combinations of states and choices before starting the value function iteration:

```julia
sol = solve(prob; rewardcall=:pre)
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

We can convert the simulation object into an array with `Array(sim)`. The order of the variables is states, choices and value (if requested).

We can also convert the simulation object into a dataframe:
```julia
df_sim = DataFrame(sim)
```

We can then use the dataframe to analyse the simulation, for example by plotting a histogram:
```julia
using Plots
histogram(df_sim.state_1, xlabel="K",legend=false)
```

or the each firm's capital stock over time:
```julia
using StatsPlots
@df df_sim plot(:period, :state_1, group=:firm, xlabel="time", ylabel="K")
```

By default, all firms start with the same state variables in the simulation. The [simulator options section](#Simulator-Options-1) contains more details about this and more sophisticated methods to intialize the simulation.

# Problem Options

## Lower integration dimension for speed

The solver must calculate expectations of future states. By default, it must integrate over all state variables, choice variables and shocks, i.e. ``s_{t+1} = \Phi(s_t, a_t, \varepsilon_{t+1})``. Many problems can be rewritten in a way such that we can integrate over fewer variables, which leads to orders of magnitude speedup.

The first step is to rewrite the problem such that the action is equal to the corresponding next period's state. In our example, the action ``a`` of the firm then is not the investment rate but simply next period's capital stock ``K'``.

```math
V(K,z) = \max_a K^\alpha e^z - i K - \frac{\gamma}{2} i^2 K + \beta \mathbb{E} V(K', z') \\
i \equiv \frac{a}{K} - (1-\delta) \\
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


Summarizing he four different types of integration dimensions:

- `:All` - next period states as a function of states, actions and shocks, i.e. ``s' = \Phi(s, a, \varepsilon')``. This formulation is the slowest and corresponds to the baseline example. In the background, the solver calls [QuanEcon's ddp.jl](https://github.com/QuantEcon/QuantEcon.jl/blob/master/src/markov/ddp.jl).
- `:Separable` - next period exogenous states as a function of all states, actions and shocks, i.e. ``s'_e = \Phi(s, a, \varepsilon')``.
- `:Separable_States` - next period exogenous states as a function of all states and shocks, i.e. ``s'_e = \Phi(s, \varepsilon')``.
- `:Separable_ExogStates` - next period exogenous states as a function of only exogenous states and shocks, i.e. ``s'_e = \Phi(s_e, \varepsilon')``.

For the solver to be fast one should use the lowest possible dimension possible. For all separable integration dimensions one must specify which state vector corresponds to the choice vector, i.e. `tChoiceVectors = (1,)`. Important: when defining the state space `tStateVectors`, the endogenous state variable(s) must come before the exogenous variables.

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

From experience, these two options can lead to orders of magnitude speedup.

The `monotonicity` keyword allows the solver to exploit the monotonicity of the choice of next period's state in the current state. Put differently, if ``s^{d*}_{t+1}`` is the optimal choice for ``s^{d*}_t``, then the optimal choice for ``s^{d'}_t \geq s^{d*}_t`` is ``s^{d'}_{t+1} \geq s^{d*}_{t+1}`` (keeping fixed the other state variables). The default is `false`.

The `concavity` keyword allows the solver to exploit the concavity of the value function in the choice of next period's state. Put differently, if ``V(s^{d'}_{t+1}, s^d_t, s^s_t) < V(s^{d*}_{t+1}, s^d_t, s^s_t)``, then also ``V(s^d_{t+1}, s^d_t, s^s_t) < V(s^{d'}_{t+1}, s^d_t, s^s_t)`` for any ``s^d_{t+1} > s^{d'}_{t+1}``. The default is `false`.

A handy way to check whether the problem fulfils these conditions is the `isapprox` function. It checks whether two different solutions are (almost) identical:

```julia
sol = solve(prob; concavity=false)
sol_conc = solve(prob; concavity=true)
isapprox(sol, sol_conc; rtol=1e-4)
```

Note: `monotonicity` and `concavity` currently only work if the integration dimension is separable.

The options `monotonicity` and `concavity` typically only work if the reward function is continuous (and monotone). In some models the reward function is monotone apart from a discontinuity at one point. For example, we could assume that the firm incurs a cost proportional to its current capital stock when it chooses a different capital stock next peried, i.e. `K' \nq K`. In this case, it is possible to still use `monotonicity` and `concavity` if we add the following option in the problem definition. We must specify that we want to try an additional value of vChoices when finding the optimal choice for a given state:

```julia
function get_additional_index(vStatesIndex) # vStatesIndex is the state variable index
	vChoice = vStatesIndex[1] # check Kprime = K
	return vChoice
end
prob = DiscreteDynamicProblem(
    tStateVectors,
    tChoiceVectors,
    reward,
    transition,
    shockdistribution,
    β;
    get_additional_index = get_additional_index)
```

This option only works for single choice variable problems at the moment.


### Prebuilding the reward matrix

The `rewardcall` keyword determines whether (part of) the reward for each state-action pair should be precomputed before the value function iteration. Prebuilding requires more memory but less computations during the VFI.

`:jit` - the reward is calculated during the VFI whenever it is required, i.e. just in time. This requires little memory but more computations during the VFI. This is the default.

`:pre` - precompute the reward for each possible combination of states and actions before the VFI. This requires more memory but less computations during the VFI.

`:pre_partial` - precompute part of the reward before the VFI that only depends on states, but not choices. From experience, this option is the fastest when combined with monotonicity and concavity. To exploit this option, we must supply the inner and outer part of the reward function when defining the problem. The outer part is the same argument in `DDP` as the standard reward function and the partial reward function enters as a keyword arguement `rewardfunc_partial`. In the solver we must then specify that the reward should be partially precomputed.

```julia
function reward(partial_reward, vStates, vChoices)
    K = vStateVars[1]
    Kprime = vChoices[1]
    i = Kprime/K - (1-δ)
    return partial_reward - γ/2*i^2 * K) - i*K
end
reward_partial(vStates) = vStates[1]^α * exp(vStates[2])

prob = DiscreteDynamicProblem(
    tStateVectors,
    tChoiceVectors,
    reward,
    transition,
    shockdistribution,
    β;
    rewardfunc_partial = reward_partial)

solve(prob; rewardcall = :pre_partial)
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
- `max_iter`: Maximum number of iterations before stopping. Default is 500.
- `disp`: whether to display number of iterations and epsilon-convergence during solving. Default is false.
- `disp_each_iter`: wait how many iterations for displaying status during solving. Default is 10.

# Simulator Options

The keyword `get_value` determines whether the simulator should also compute the value at the beginning of each time period. This takes more time. The default is `get_value = false`.

## Initialization

The default setting of how the initial state variables are drawn at t=0 is the burn-in method. Here, each firm starts with states that are the median values of the state space. The idea is to simulate the model for more periods than necessary and then only analyse the simulated values after 60 periods or so on, depending on the persistence in the model.

There is also an option for more sophisticated initialization, where the initial endogenous state variables are chosen by the firm subject to a loss function and the exogenous ones are predetermined or random. In the example, we assume the initial state variables are determined as follows:

```math
V_0(z_0) = \max_{K_0} V(K_0,z_0) - (1 + (1 + C0) * K_0 \\
\text{where } z_0 = \sqrt{\frac{\sigma^2}{1-\rho^2}} ε_0 \\
ε_0 \sim \mathcal{N}(0,1)
```

The initial productivity $z_0$ is drawn from a normal distribution. Depending on $z_0$, the firm chooses it's initial capital stock $K_0$ to maximize it's continuation value, accounting for the fact that there is a deadweight cost $C_0$ to acquire $K_0$. If $C_0 > 0$, then $K_0$ will be optimally chosen below it's steady state value.

We must code the following as an additional input when constucting `DiscreteDynamicProblem`:

```julia
initprob(value, vChoices) = value - (1 + (1-β)/β + C0) * vChoices[1]
init(vShocks) = vShocks[1] * sqrt(σ^2 / (1-ρ^2)) # = z0
prob = DiscreteDynamicProblem(...;
    initializationproblem = initprob,
    initializefunc = init,
    shockdist_initial = Normal(),
    tChoiceVectorsZero = (1,)) # which state variable does the firm choose at t=0
```

The `solve` function then also finds the optimal policy at t=0 when the problem is defined this way. We can access the solutions via:

```julia
sol = solve(prob)
value0(sol)[5] # initial value if z0 = fifth grid point
value0(sol)(0.5) # initial value if z0 = 0.5
policy0(sol)[5] # initial policy K_0 if z0 = fifth grid point
policy0(sol)[5, 3] # initial policy K_0 if z0 = 0.5
```

Even if the initializationproblem is specified this way, there is still an option to resort to the burn-in method. The keyword argument `initialize_exact` controls whether the simulator should use the sophisticated method or not. The default is `true`.

## Fixed shock draws

It is also possible to simulate a model solution using an identical draw of shocks:

```julia
shocks = drawshocks(prob, nPeriods = 60, nFirms = 100)
sim = simulate(sol, shocks)
```

This can be useful for doing comparative statics on the model parameters, and wants to be sure that the shocks do not change. Note: the draw of shocks refers to the supplied shock distribution in the problem defintion. If the shock distribution is parametrized, for example by its variance, then one should not do comparative statics on those parameters.
