# DiscreteDynamicProgramming.jl Documentation

This package solves and simulates discrete dynamic choice models with value function iteration fast. By taking advantage of properties of the problem, a wide array of problems can be solved very fast without requiring parallelization.

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

We code a neoclassical capital investment model with convex adjustment costs as an example. It's code is as follows:

```math
V(K,z) = \max_a K^\alpha e^z - (a - (1-\delta) K) + \beta \mathbb{E} V(K', z') \\
\text{where } K' = a \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
```

To code this model:

```julia
β = 0.9 # discount factor
α = 0.67 # returns to scale parameter 
ρ = 0.6 # persistence of producitivity
σ = 0.3 # volatility of productivity
δ = 0.15 # deprectiation rate of capital
γ = 2.0 # convex adjustment cost parameter

function rewardfunc(vStateVars, vChoices)
    (K, z) = vStateVars
    Kprime = vChoices[1]
    capx = Kprime - (1-δ)*K
    return K^α * exp(z) + γ/2*(capx/K- δ)^2 * K
end

function transfunc(vStateVars, vChoices, vShocks)
    z = vStateVars[2]
    Kprime = vChoices[1]
    zprime  = ρ*z + σ * vShocks[1];
    return  Kprime, zprime
end

shockdistribution = Normal() # standard normal

# create state variable vectors
nz = 5
stdz = sqrt(σ^2/(1-ρ^2))
minz = -3*stdz
maxz =  3*stdz
vz = collect(LinRange(minz, maxz, nz))

nK = 100
logK_ss = log(α * exp(0.)/ ((1-β)/β*(1+γ*δ) + δ)) ^ (1/(1-α))
log_K_min = 0.1 * logK_ss
log_K_max = 2 * logK_ss
vK   = exp.(collect(LinRange(log_K_min, log_K_max, nK)))

tStateVectors = (vK, vz)
tChoiceVectors = (vK,)

prob = DiscreteDynamicProblem(
            β,
            rewardfunc,
            transfunc,
            shockdistribution,
            tStateVectors,
            tChoiceVectors,
            )
```

## Step 2: solve problem

Solving it with the following, which provides a mesh of the value function and policy function.

```julia
sol = solve(prob)
```

There are different options available which make solving the model much faster.

#### Options

`intdim`: specify the integration dimension of the problem. By default, we integrate over all dimensions of the state variables. Setermines over how many variables we integrate when calculating expectations. There are three possible options.

* `:SA`: integrate across each possible combination of state variables and choice variables. Assumes that shocks, current states and current choices affect transition probabilities of the stochastic variables, i.e. ``s_{t+1} = \Phi(s_t, a_t, \varepsilon_{t+1})``, where ``\varepsilon_{t+1} \sim \mathcal{N}(0,1)``. This is the most general but slow (is implemented by the QuantEcon library). It's good for testing with a small state-space.

The other two options require that we can separate the state variables into two groups, exogenous and endogenous variables. The transition of the stochastic state variables ``s^s`` only depends on shocks and states but not actions, i.e.  ``s^s_{t+1} = \Phi(s_t, \varepsilon_{t+1})``. The transition of deterministic state variables ``s^d`` only depends on states and actions but not shocks, i.e.  ``s^d_{t+1} = \Phi(s_t, a_t)``. To operationalize this, we substitute the transition of deterministic state variables into the reward, i.e. ``r(s_t, a_t) = r(s_t, \Phi^{-1}(s_t, s^d_{t+1})) ``. In that sense, our choice is in terms of choosing the future deterministic state directly.

* `:intermediate`: only states, but not current choices affect the transition of the stochastic state variables. For example, in a version of the R&D model, the R&D stock inherited from last period affects the transition of productivity, but not current R&D expenditures. Change the transition function like so:
```julia
function transfunc(vStateVars, vShocks)
    z = vStateVars[2]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
```
and mention in the options that `intdim = :intermediate`.


* `:separable`: is when neither deterministic state variables nor choices affect the transition of the stochastic state variables, i.e. ``s^s_{t+1} = \Phi(s^s_t, \varepsilon_{t+1})``. For example in the Neoclassical model, the capital stock does not affect future productivity. Change the transition function like so:
```julia
function transfunc(vExogStateVars, vShocks)
    z = vStateVars[1]
    zprime  = ρ*z + σ * vShocks[1];
    return  zprime
end
```
and mention in the options that `intdim = :separable`.


The options `monotonicity` and `concavity` exploit properties of the value function to speed up the VFI. They don't work with `:intdim = :SA`.

`monotonicity`: exploit the monotonicity of the choice of next period's state in the current state. Put differently, if ``s^{d*}_{t+1}`` is the optimal choice for ``s^{d*}_t``, then the optimal choice for ``s^{d'}_t \geq s^{d*}_t`` is ``s^{d'}_{t+1} \geq s^{d*}_{t+1}`` (keeping fixed the other state variables).

`concavity`: exploit the concavity of the value function in the choice of next period's state. Put differently, if ``V(s^{d'}_{t+1}, s^d_t, s^s_t) < V(s^{d*}_{t+1}, s^d_t, s^s_t)``, then also ``V(s^d_{t+1}, s^d_t, s^s_t) < V(s^{d'}_{t+1}, s^d_t, s^s_t)`` for any ``s^d_{t+1} > s^{d'}_{t+1}``.

`rewardmat`: determines whether (part of) the reward for each state-action pair should be pre-built before the value function iteration. Prebuilding it requires more memory but less computations during the VFI. Some parts of it may never be used with concavity and monotonicity.

* `:nobuild` means that the reward is calculated during the VFI. This requires little memory but more computations during the VFI.
* `:prebuild_partial` means that part of the reward function is prebuilt, the one that does not depend on choices, but only states. From experience, this is the fastest.
* `:prebuild` means prebuild the entire reward matrix for each possible combination of states and actions. This requires more memory but less computations during the VFI.

The solver calculates the transition matrix for you using Gaussian Quadrature. The user must supply the transition function. The shocks are assumed to be standard random normal. Optionally, the transition matrix can also be provided manually,

```julia
using QuantEcon: tauchen
mTransition = tauchen(prob.params.nz, prob.params.ρ, prob.params.σ)
sol = solve(prob, mTransition = mTransition)
```

The option `disp = true` determines whether the iterations should be printed in the REPL.


## Step 3: simulate the model

To simulate a panel:

```julia
nPeriods = 120
nFirms = 1000
sim = simulate(prob, sol, nPeriods, nFirms)
```

If `nPeriods` and `nFirms` are contained in the `prob.params` structure, then we can simulate directly with `simulate(prob, sol)` (if the information is provided, then it is overridden). It is also possible to simulate a model or variations thereof (e.g. for different parameters) using the same draw of shocks:

```julia
nPeriods = 120
nFirms = 1000
shocks = drawshocks(prob, nPeriods, nFirms)
sim = simulate(prob, sol, shocks)
```

# notes

## Solver assumptions

- state variables are always forced to stay within bounds
- for `intermediate` and `separable`, can choose directly the state variable
- order of state variables is first the endogenous state variable, then exogenous

## Reserved parameters

- discount factor beta
- solver options
