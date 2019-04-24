using DiscreteDynamicProgramming

function warn_time(t, border, preprint)
	t < border || @warn preprint * " took longer than usual: $t vs $border seconds"
end
function warn_memory(bytes, border, preprint)
	bytes < border || @warn preprint * " allocated more memory than usual: " *
		Base.format_bytes(bytes) * " vs. " * Base.format_bytes(border)
end


# Single choice variable
model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5, :F => 0., :τ => 0.3, :β => 0.9)
p_neoclassical_All = createmodel(model; dipar..., intdim = :All)
p_neoclassical_Separable = createmodel(model; dipar..., intdim = :Separable)
p_neoclassical_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_neoclassical_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardcall=>:pre)

sol_neoclassical_All = solve(p_neoclassical_All; diopt...)
val, t, bytes = @timed solve(p_neoclassical_All; diopt...)
warn_time(t, 0.12, "Neo All: Solver")
warn_memory(bytes, 75555555, "Neo All: Solver")

sol_neoclassical_Separable = solve(p_neoclassical_Separable; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable; diopt...)
warn_time(t, 0.07, "Neo Separable: Solver")
warn_memory(bytes, 19000000, "Neo Separable: Solver")

sol_neoclassical_Separable_States = solve(p_neoclassical_Separable_States; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable_States; diopt...)
warn_time(t, 0.01, "Neo Separable_States: Solver")
warn_memory(bytes, 650000, "Neo Separable_States: Solver")

sol_neoclassical_Separable_ExogStates = solve(p_neoclassical_Separable_ExogStates; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable_ExogStates; diopt...)
warn_time(t, 0.01, "Neo Separable_ExogStates: Solver")
warn_memory(bytes, 210000, "Neo Separable_ExogStates: Solver")

# SIMULATE
shocks = drawshocks(p_neoclassical_All, nFirms=10^4, nPeriods=50)

simulate(sol_neoclassical_All, shocks)
val, t, bytes = @timed simulate(sol_neoclassical_All, shocks)
warn_time(t, 0.25, "Neo All: Simulator")
warn_memory(bytes, 221000000, "Neo All: Simulator")

simulate(sol_neoclassical_Separable, shocks)
val, t, bytes = @timed simulate(sol_neoclassical_Separable, shocks)
warn_time(t, 0.25, "Neo Separable: Simulator")
warn_memory(bytes, 260000000, "Neo Separable: Simulator")

simulate(sol_neoclassical_Separable_States, shocks)
val, t, bytes = @timed simulate(sol_neoclassical_Separable_States, shocks)
warn_time(t, 0.24, "Neo Separable_States: Simulator")
warn_memory(bytes, 260000000, "Neo Separable_States: Simulator")

simulate(sol_neoclassical_Separable_ExogStates, shocks)
simulate(sol_neoclassical_Separable_ExogStates, shocks)
val, t, bytes = @timed simulate(sol_neoclassical_Separable_ExogStates, shocks)
warn_time(t, 0.23, "Neo Separable_ExogStates: Simulator")
warn_memory(bytes, 190000000, "Neo Separable_ExogStates: Simulator")


# Two choice variables
model = :Intangible
dipar = Dict(:nK => 10, :nN => 10, :nz => 5, :β => 0.9)
p_int_All = createmodel(model; dipar..., intdim = :All)
p_int_Sep = createmodel(model; dipar..., intdim = :Separable)
p_int_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_int_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardcall=>:pre)

sol_int_All = solve(p_int_All; diopt...)
val, t, bytes = @timed solve(p_int_All; diopt...)
warn_time(t, 2.1, "Intangible All: Solver")
warn_memory(bytes, 860000000, "Intangible All: Solver")

sol_int_Sep = solve(p_int_Sep; diopt...)
val, t, bytes = @timed solve(p_int_Sep; diopt...)
warn_time(t, 1., "Intangible Separable: Solver")
warn_memory(bytes, 800000000, "Intangible Separable: Solver")

sol_int_Sep_States = solve(p_int_Sep_States; diopt...)
val, t, bytes = @timed solve(p_int_Sep_States; diopt...)
warn_time(t, 0.1, "Intangible Separable_States: Solver")
warn_memory(bytes, 2000000, "Intangible Separable_States: Solver")

sol_int_Sep_ExogStates = solve(p_int_Sep_ExogStates; diopt...)
val, t, bytes = @timed solve(p_int_Sep_ExogStates; diopt...)
warn_time(t, 0.1, "Intangible Separable_ExogStates: Solver")
warn_memory(bytes, 552000, "Intangible Separable_ExogStates: Solver")

# SIMULATE
shocks = drawshocks(p_int_All, nFirms=10^3, nPeriods=50)

simulate(sol_int_All, shocks)
val, t, bytes = @timed simulate(sol_int_All, shocks)
warn_time(t, 0.12, "Int All: Simulator")
warn_memory(bytes, 62000000, "Int All: Simulator")

simulate(sol_int_Sep, shocks)
val, t, bytes = @timed simulate(sol_int_Sep, shocks)
warn_time(t, 0.07, "Int Separable: Simulator")
warn_memory(bytes, 48000000, "Int Separable: Simulator")

simulate(sol_int_Sep_States, shocks)
val, t, bytes = @timed simulate(sol_int_Sep_States, shocks)
warn_time(t, 0.07, "Int Separable_States: Simulator")
warn_memory(bytes, 48000000, "Int Separable_States: Simulator")

simulate(sol_int_Sep_ExogStates, shocks)
val, t, bytes = @timed simulate(sol_int_Sep_ExogStates, shocks)
warn_time(t, 0.07, "Int Separable_ExogStates: Simulator")
warn_memory(bytes, 48000000, "Int Separable_ExogStates: Simulator")
