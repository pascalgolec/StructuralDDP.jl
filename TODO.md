- [c] allow discontinuities in the reward function with `monotonicity` and `concavity` for two choice variables
    - no speed improvement for `monotonicity`, only minor for `concavity`
- [x] improve indexing inside VFI: count up the index in a smarter way, then don't have to calculate from scratch every time
- More example models (have two testmodels in test/models for now)
- Benchmark with VFI with GPU ([here](https://discourse.julialang.org/t/value-function-iteration-on-gpu/13774/6) and [here](https://juliacon.org/2018/talks_workshops/106/)) and [others](https://github.com/JuliaStochOpt/StochDynamicProgramming.jl)
- check [Event Handling and Callback Functions](http://docs.juliadiffeq.org/latest/features/callback_functions.html) in DifferentialEquations.jl
- support more than two choice variables when the integration dimension is `separable`
- plotrecipe for `DDPSolution`
- allow parametric types in problem definition of state and choice vectors, i.e. `AbstractVector{T}` instead of `Vector{Float64}`
- allow `LabelledArrays.jl` in problem definition
- support firm exit
- add @inbounds
- fixed shock draws should be possible if [fix the local RNG](https://github.com/JuliaStats/Distributions.jl/issues/436)
- :firm in DataFrame(sim) should be an CategoricalArray for easy plotting
- document plotting with VegaLite.jl
```julia
f_sim.firm = CategoricalArray(df_sim.firm)
plt = df_sim|> @vlplot(
    :bar,
    x={:state_1, bin=true, title="capital stock"},
    y="count()")
plt = df_sim|> @vlplot(
    :line,
    x={:period, title="years"},
    y={:state_1, title="capital stock"},
    color=:firm)
```
