# DiscreteDynamicProgramming.jl Documentation

This package solves and simulates discrete dynamic choice models with value function iteration. It has a simple and transparent syntax for defining a dynamic optimization problem and is optimized for speed.

It includes options to take advantage of properties of the problem which can result in orders of a magnitude speedup. Solving the models fast without requiring parallelization is useful if one ultimately wants to simulate the model for different parameter values, as then one can parallelize across that dimension.

The solver is fast for problems where some of next periods state variables can be chosen directly.

The general workflow is to define a problem, solve the problem, and then analyze the solution. The full code for solving and simulating the example is:

```julia
using DiscreteDynamicProgramming
prob = NeoClassical(nK=150, nz=15, ρ=0.5, σ=0.3, γ=2.)
sol = solve(prob, intdim=:separable,
    monotonicity=true, concavity=true,
    rewardmat=:prebuild_partial)
sim = simulate(prob, sol)
```

where the pieces are described below.

## Step 1: define problem

Each discrete dynamic maximization problem has an objective function as follows

```math
\mathbb{E} \sum_{t=0}^{\infty} \beta^t r(s_t, a_t)
```

where

- ``r(s_t, a_t)`` is the current reward as a function of the state variables ``s_t`` and choice variable(s) ``a_t``
- transition equation(s) for state variables ``s_{t+1}`` depending on ``(s_t, a_t)`` and i.i.d shocks ``\varepsilon_{t+1}``.
- a discount factor ``β``

For a more formal definiton, see the QuantEcon lectures [here](https://lectures.quantecon.org/jl/discrete_dp.html#Discrete-DPs). We can write the problem in recursive form:

```math
V(s) = \max_a r(s, a) + \beta \mathbb{E} V(s') \\
\text{where } s' = \Phi(s, a, \varepsilon')
```

where ``s'`` corresponds to next period's state ``s_{t+1}``.

Note: also mention closed action and state space. Look in notation of Acemoglu book.

We will solve a capital investment model with convex adjustment costs as an example. It is defined as follows:

```math
V(K,z) = \max_I K^\alpha e^z - I - \frac{\gamma}{2} \frac{I^2}{K} + \beta \mathbb{E} V(K', z') \\
\text{where } K' = (1-\delta)K + I \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

where ``K`` is current and ``K'`` next period's capital, ``I`` is capital investment. Capital depreciates at a rate ``\delta``. Investment ``I`` adds to the capital stock which becomes productive next period. The parameter ``\gamma`` determines the adjustment costs for investment, over and above the cost of capital ``I``. The productivity of capital is determined by ``z``, which follows an AR(1) process with autocorrelation ``\rho`` and volatility ``\sigma``. There are decreasing returns to scale for ``\alpha < 1``, so there is an optimal scale of the firm depending on ``z``. So the problem of the firm is to choose the optimal scale for next period, taking into account that ``z'`` is unknown and large adjustments are disproportianally costly.

We define the problem as follows:

```julia
β = 0.9 # discount factor
α = 0.67 # returns to scale parameter
ρ = 0.6 # persistence of producitivity
σ = 0.3 # volatility of productivity
δ = 0.15 # deprectiation rate of capital
γ = 2.0 # convex adjustment cost parameter


# create state variable vectors
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
            β,
            rewardfunc,
            transfunc,
            shockdistribution,
            tStateVectors,
            tChoiceVectors,
            )
```

In some problems one may only have one shock or one choice variable. In this case, one can write the reward- and transition functions more concisely in terms as numbers as inputs as opposed to vectors (it's actually slighlty faster):

```julia
function reward(vStateVars, Kprime)
    (K, z) = vStateVars
    capx = Kprime - (1-δ)*K
    return K^α * exp(z) + γ/2*(capx/K- δ)^2 * K
end
function transition(vExogState, shock)
    z = vExogState[1]
    zprime  = ρ*z + σ * shock;
    return inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
end
```

## Step 2: solve problem

Solving it with the following, which provides a mesh of the value function and policy function.

```julia
sol = solve(prob)
```

The solver can be controlled using different options which are discribed in the [Solver Options section](#Solver-Options-1). For example, we can display the value function iterations with `disp`:

```julia
sol = solve(prob; disp=true)
```

## Step 3: simulate the model

We can simulate a panel:

```julia
sim = simulate(prob, sol, nPeriods = 120, nFirms = 10000)
```

It is also possible to simulate a model or variations thereof using the identical draw of shocks:

```julia
shocks = drawshocks(prob, nPeriods = 120, nFirms = 1000)
sim = simulate(prob, sol, shocks)
```

Note: the draw of shocks refers to the supplied shock distribution in the problem defintion. If the shock distribution is parametrized, for example by its variance, then one should not do comparative statics on those parameters.

### Initialization

The default setting how the initial state variables are drawn at t=0 is the burn-in method. Here, each firm starts with states that are the median values of the state space. The idea is to simulate the model for more periods than necessary and then only analyse the simulated values after 60 periods or so on, depending on the model.

There is also an option for more sophisticated initialization, where the initial endogenous state variables are chosen by the firm subject to a loss function and the exogenous ones are predetermined or random. For the neoclassical model, this entails coding the following functions as an additional input into the problem:

```julia
initprob(value::Float64, vChoices) = value - (1 + (1-β)/β + C0) * vChoices[1]
function init(vShocks)
    z0 = vShocks[1] * sqrt(σ^2 / (1-ρ^2))
    return inbounds(z0, tStateVectors[2][1], tStateVectors[2][end])
end
createDiscreteDynamicProblem(<other variables>,
    initializationproblem = initprob,
    initializefunc = init,
    tChoiceVectorsZero = (1,))
```

If we have a problem where our choice variables are not exacly equal to some of the state variables, then we need to specify additionally which state variables the firm chooses at t=0. This is achieved with the option `tChoiceVectorsZero`, i.e. `tChoiceVectorsZero = (1,)` for the neoclassical model. If `tChoiceVectorsZero` is not specified, then the solver/simulator uses indices of the choice variables for the dynamic optimization.

# Solver Options

## For speed

There are different options available which make solving the model much faster.

### Lower integration dimension

The solver is fast if we slightly rewrite the problem, where the action ``a`` of the firm is not investment but simply next period's capital stock ``K'``.

```math
V(K,z) = \max_a K^\alpha e^z - (a - (1-\delta) K) - \frac{\gamma}{2} \frac{(a - (1-\delta) K)^2}{K} + \beta \mathbb{E} V(K', z') \\
\text{where } K' = a \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

The transition function generally is a function of the state variables, choices and shocks. In the captal investment example however where the action is equal to next period's state, we can sort the state variables into endogenous and exogenous state variables, ``K`` and ``z``. What's more, neither ``K`` nor the choice affects the transition of ``z``. Thus we can define the transition function in a compact manner as follows:

By default, we integrate over all state variables, choice variables and shocks to get expectations of next periods states, i.e. ``s_{t+1} = \Phi(s_t, a_t, \varepsilon_{t+1})``. Many problems however don't require this.

In the Neoclassical model, we can directly choose next period's `K` and don't need to integrate along that dimension.
Formally, this is the case when the choice variable is equal to one of the state variables. Note: sometimes this requires rewriting the model by redefing the action space.

In that case, we can provide that information when defining the problem by specifying which state variable is the choice variable instead of defining directly the choice vectors, i.e. `tChoiceVectors = (1,)` which means that the first choice variable is equal to the first state variable. Important: when defining the state space `tStateVectors`, the endogenous state variable(s) must come first. We also need to rewrite the transition function to only return `z` (note that in principle `zprime` could depend on the choice of `Kprime`:
```julia
function transfunc(vStateVars, vChoices, vShocks)
    z = vStateVars[2]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
```
We also need to mention in the problem definition which choice variable corresponds to the index of which state variable.

Formally, we require that we can separate the state variables into two groups, endogenous and exogenous state variables, which correspond to `K` and `z` in the example.
The transition of the endogenous state variables ``s^d`` only depends on states and actions but not shocks, i.e.  ``s^d_{t+1} = \Phi(s_t, a_t)``.
The transition of the stochastic state variables ``s^s`` only depends on shocks and states but not actions, i.e.  ``s^s_{t+1} = \Phi(s_t, \varepsilon_{t+1})``. Note: sometimes this requires slightly rewriting the model by redefing the action space.

If this separation is possible, then then change the transition function like so:
```julia
function transfunc(vStateVars, vShocks)
    z = vStateVars[2]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
```

and mention it in the model options, `intdim = :intermediate`.

The solver is even faster if additionally only the shocks and exogenous state variables, but not endgenous state variables, determine next period's exogenous state variables, i.e. ``s^s_{t+1} = \Phi(s^s_t, \varepsilon_{t+1})``. In that case, change the transition function like so (note the change in the index of `z`):
```julia
function transfunc(vExogStateVars, vShocks)
    z = vExogStateVars[1]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
```
and mention in the model options, `intdim = :separable`.

bEndogStateVars?

### Monotonicity and concavity

The options `monotonicity` and `concavity` exploit properties of the value function to speed up the VFI. They don't work with `:intdim = :SA`.

`monotonicity`: exploit the monotonicity of the choice of next period's state in the current state. Put differently, if ``s^{d*}_{t+1}`` is the optimal choice for ``s^{d*}_t``, then the optimal choice for ``s^{d'}_t \geq s^{d*}_t`` is ``s^{d'}_{t+1} \geq s^{d*}_{t+1}`` (keeping fixed the other state variables).

`concavity`: exploit the concavity of the value function in the choice of next period's state. Put differently, if ``V(s^{d'}_{t+1}, s^d_t, s^s_t) < V(s^{d*}_{t+1}, s^d_t, s^s_t)``, then also ``V(s^d_{t+1}, s^d_t, s^s_t) < V(s^{d'}_{t+1}, s^d_t, s^s_t)`` for any ``s^d_{t+1} > s^{d'}_{t+1}``.


### Prebuilding the reward matrix

`rewardmat`: determines whether (part of) the reward for each state-action pair should be pre-built before the value function iteration. Prebuilding it requires more memory but less computations during the VFI. Some parts of it may never be used with concavity and monotonicity.

* `:nobuild` means that the reward is calculated during the VFI. This requires little memory but more computations during the VFI.
* `:prebuild_partial` means that part of the reward function is prebuilt, the one that does not depend on choices, but only states. From experience, this is the fastest.
* `:prebuild` means prebuild the entire reward matrix for each possible combination of states and actions. This requires more memory but less computations during the VFI.

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

# notes

## Solver assumptions

- state variables are always forced to stay within bounds
- for `intermediate` and `separable`, can choose directly the state variable
- order of state variables is first the endogenous state variable, then exogenous
