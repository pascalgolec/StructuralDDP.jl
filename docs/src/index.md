# DiscreteDynamicProgramming

This package solves and simulates discrete dynamic choice models with value function iteration fast. By taking advantage of properties of the problem, a wide array of problems can be solved very fast without requiring parallelization. Inline latex: ``\sqrt{1 + x + x^2 + }``

Here's an equation:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

The general workflow is to define a problem, solve the problem, and then analyze the solution. The full code for solving an example is:

```julia
using DiscreteDynamicProgramming
prob = createmodel(:NeoClassicalSimple, nK=150, nz=15, ρ=0.5, σ=0.3, γ=2.)
sol = solve(prob, :intdim=>:separable, :monotonicity=>true, :concavity=>true, :rewardmat=>:prebuild_partial)
```

## Problem definition

Each discrete dynamic maximization problem has an objective function as follows

```math
\mathbb{E} \sum_{t=0}^{\infty} \beta^t r(s_t, a_t)
```

where

- ``r(s_t, a_t)`` is the current reward as a function of the state variables ``s_t`` and choice(s) ``a_t``
- transition equation(s) for state variables ``s_{t+1}`` depending on ``(s_t, a_t)``
- a discount factor ``β``

For a more formal definiton, see the QuantEcon lectures [here](https://lectures.quantecon.org/jl/discrete_dp.html#Discrete-DPs).

The package includes a neoclassical capital investment model with convex adjustment costs as an example. We define the problem by creating an instance of the model with specific parameters,

```julia
prob = createmodel(:NeoClassicalSimple, ρ=0.7, σ=0.2, β=0.9)
```

To solve the model,

```julia
sol = solve(prob, :intdim=>:separable, :monotonicity=>true, :concavity=>true, :rewardmat=>:prebuild_partial)
```

### Preparation

Need to prepare some functions if want to estimate your own model.

- definition of the Type and the constructor, in the form of `createmodel(Type)`
- To solve, need
	- `rewardfunc`, that gives the flow reward as a function of the state variables and choice variables
		- if you use nested functions you can also use them when calculating the accounting variables from the simulated model
	- `transfunc`, a function of states, choices and shocks, to get the evolution of the state variables
		- nest them if possible?
	- `initializationproblem`, needed to construct matrices of indirect utility and initial choice depending on exogenous variables. Use it later to initialize simulation
	- for the `separate` solver, additionally need
		- `outputfunc`, to calculate output from the current state in the form of a matrix

- To simulate
	- `initialize`, how the simulation should be initalized
- To calculate accounting variables
	- `accountingvars`, takes array of simulated state variables and constructs a dataframe with accounting variables
	- `dividends`, specify if dividends are unequal to profits
	- `adjustcosts`, specify if have other than convex and fixed adjustment costs
	- `oibdp`, specify if not calculated as grossprofits - adjustcosts
	- `grossprofits`, output as function of state vars

### Solver

There are different solvers for the models. When creating an instance of a model, specify which solver you'd like to use.

- Support up to two choice variables for now.
- Need that choice vector equal to state vector I think?

```julia
ptest = createmodel(NeoClassicalSimple, intdim = :SA, monotonicity = true, concavity = false)
```
#### Options

`rewardmat`: determines whether (part of) the reward matrix should be pre-built before the VFI. Prebuilding it can require a lot of memory, and some parts of it may never be used with concavity and monotonicity.

* `:nobuild` means that the reward function is always recalculated
* `:prebuild_partial` means that part of the reward function is prebuilt, the one that does not depend on choices, but only states.
* `:prebuild` means prebuild the entire reward matrix for each possible combination of states and actions. (This is what is used in the State-Action formulation.)

`intdim`: determines over how many variables we integrate.

* The most general but slowest is the State-action (`:SA`) pair computation which is implemented by the QuantEcon library. It's good for testing with a small state-space. Current states and current choices are allowed to affect transition probabilities of the stochastic variables.

The following formulations assume that at least one of the next period states can exactly be chosen without any randomness, for example the capital equation.

* The second formulation (`:intermediate`) is where only states, but not current choices affect the transition of the stochastic state variables. For example, in a version of the R&D model, the R&D stock inherited from last period affects the transition of productivity, but not current R&D expenditures.
* A third formulation (`:separable`) is when the endogenous state variables don't affect the transition of the stochastic state variables. For example in the Neoclassical model, the capital stock does not affect future productivity.
	* This is faster because we need to integrate over less variables

`monotonicity`: exploit the monotonicity of the policy function of the first state variables (which is also a choice variable). During the VFI, as we iterate across optimal candidates of the first state variable, we start with the last optimal candidate and not the first one.

`concavity`: exploit the concavity of the value function in the first state variable (which is also a choice variable). During the VFI, as we iterate across optimal candidates of the first state variable, we break when the value decreases, as from then on all higher candidates' value will be lower.

`monotonicity` and `concavity` don't work with `:SA`.


# notes

## Solver assumptions

- state variables are always forced to stay within bounds
- for `intermediate` and `separable`, can choose directly the state variable
- order of state variables is first the endogenous state variable, then exogenous

## Reserved parameters

- discount factor beta
- solver options
