# This script solves two capital adjustment models in two different ways. Adding
# no structure to the problem versus adding as much structure as possible.
# 	- no structure: state-action formulation
# 	- structure: add as much information as possible
# todo: benchmark compared to POMDPs.jl instead of my own package
# https://stackoverflow.com/questions/52885318/porting-an-example-from-quantecon-jl-to-pomdps-jl

using StructuralDDP
using DataFrames, Gadfly
using InteractiveUtils, BenchmarkTools

# include("../test/models/CapitalAdjustModel.jl")
# include("../test/models/CapitalAdjustModel2.jl")
# createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

par_1 = Dict(:nK => 150, :nz => 5, :γ => 2., :β=>0.9)
p_1_sa = createmodel(:CapitalAdjustModel; par_1..., intdim=:All)
p_1_str = createmodel(:CapitalAdjustModel; par_1..., intdim=:Separable_ExogStates)

# precompile
solve(p_1_sa; concavity = false, monotonicity = false, rewardcall=:pre)
solve(p_1_str, concavity = true, monotonicity = true,
	rewardcall=:pre_partial)

# measure
t_1_sa = @timed solve(p_1_sa)
t_1_str = @timed solve(p_1_str, concavity = true, monotonicity = true,
	rewardcall=:pre_partial)

par_2 = Dict(:nK => 15, :nN=>15, :nz=>3, :β=>0.9)
p_2_sa = createmodel(:CapitalAdjustModel2; par_2..., intdim=:All)
p_2_str = createmodel(:CapitalAdjustModel2; par_2..., intdim=:Separable_ExogStates)

solve(p_2_sa, concavity = false, monotonicity = false, rewardcall=:pre)
solve(p_2_str, concavity = true, monotonicity = true, rewardcall=:pre_partial)

t_2_sa = @timed solve(p_2_sa,
	concavity = false, monotonicity = false, rewardcall=:pre)
t_2_str = @timed solve(p_2_str,
	concavity = true, monotonicity = true, rewardcall=:pre_partial)

test1 = "one choice variable"
test2 = "two choice variables"
method1 = "standard"
method2 = "structured"
df = DataFrame(test = test1, method = method1, time = t_1_sa[2])
	push!(df, [test1, method2, t_1_str[2]])
	push!(df, [test2, method1, t_2_sa[2]])
	push!(df, [test2, method2, t_2_str[2]])


plt = Gadfly.plot(df,
    x = :test,
    y = :time,
	color = :method,
    Scale.y_log10,
    Guide.ylabel("seconds (orders of magnitude)"),
    Guide.xlabel(nothing),
    # Coord.Cartesian(xmin=1,xmax=13.3,ymin=-0.5,ymax=4.2),
    Theme(
        guide_title_position = :left,
        colorkey_swatch_shape = :circle,
        minor_label_font = "Georgia",
        major_label_font = "Georgia"),
	)

draw(SVGJS("compare.svg"), plt)
