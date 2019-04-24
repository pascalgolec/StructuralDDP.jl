using DiscreteDynamicProgramming, DataFrames, Gadfly
using InteractiveUtils, BenchmarkTools

include("../test/models/Neoclassical_user.jl")
include("../test/models/Intangible_user.jl")
createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

model_par = Dict(:nK => 300, :nz => 15, :γ => 2., :β=>0.95)
p_neo = createmodel(:NeoClassicalSimple; model_par..., intdim=:Separable_ExogStates)

# precompile
solve(p_neo; concavity = false, monotonicity = false, rewardcall=:pre)
solve(p_neo, concavity = true, monotonicity = true,
	rewardcall=:pre_partial)

# measure
t_neo = @timed solve(p_neo; concavity = false, monotonicity = false, rewardcall=:pre)
t_neo_fast = @timed solve(p_neo, concavity = true, monotonicity = true,
	rewardcall=:pre_partial)

model_par = Dict(:nK => 50, :nN=>25, :nz=>5, :γ=>2., :β=>0.95)
p_intan = createmodel(:Intangible; model_par..., intdim=:Separable_ExogStates)

solve(p_intan, concavity = false, monotonicity = false, rewardcall=:pre)
solve(p_intan, concavity = true, monotonicity = true, rewardcall=:pre_partial)

t_intan = @timed solve(p_intan,
	concavity = false, monotonicity = false, rewardcall=:pre)
t_intan_fast = @timed solve(p_intan,
	concavity = true, monotonicity = true, rewardcall=:pre_partial)

test1 = "one choice var"
test2 = "two choice var"
method1 = "standard"
method2 = "accelerated"
df = DataFrame(test = test1, method = method1, time = t_neo[2])
	push!(df, [test1, method2, t_neo_fast[2]])
	push!(df, [test2, method1, t_intan[2]])
	push!(df, [test2, method2, t_intan_fast[2]])


plt = Gadfly.plot(df,
    x = :test,
    y = :time,
	color = :method,
    Scale.y_log10,
    Guide.ylabel("seconds"),
    Guide.xlabel(nothing),
    # Coord.Cartesian(xmin=1,xmax=13.3,ymin=-0.5,ymax=4.2),
    Theme(
        guide_title_position = :left,
        colorkey_swatch_shape = :circle,
        minor_label_font = "Georgia",
        major_label_font = "Georgia"),
	)

draw(SVGJS("compare.svg"), plt)
