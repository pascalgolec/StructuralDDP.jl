
using DiscreteDynamicProgramming, DataFrames, Gadfly

p_neoclassical = createmodel(:NeoClassicalSimple; nK = 300, nz = 15, γ = 2., F = 0., β=0.95)
t_neo = @timed solve(p_neoclassical, intdim = :separable,
	concavity = false, monotonicity = false,
	rewardmat=:prebuild)
t_neo_fast = @timed solve(p_neoclassical, intdim = :separable,
	concavity = true, monotonicity = true,
	rewardmat=:prebuild_partial)

p_intan = createmodel(:Intangible; nK = 50, nN = 25, nz=5, β=0.95)
t_intan = @timed solve(p_intan, intdim = :separable,
	concavity = false, monotonicity = false, rewardmat=:prebuild)
t_intan_fast = @timed solve(p_intan, intdim = :separable,
	concavity = true, monotonicity = true, rewardmat=:prebuild_partial)

test1 = "neoclassical model"
test2 = "intangibles model"
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
        major_label_font = "Georgia"
    ),
	)

draw(SVGJS("compare.svg"), plt)
