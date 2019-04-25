using StructuralDDP

# load models
include("models/CapitalAdjustModel.jl")
include("models/CapitalAdjustModel2.jl")
createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

function warn_time(t, border, preprint)
	t < border || @warn preprint * " took longer than usual: $t vs $border seconds"
end
function warn_memory(bytes, border, preprint)
	bytes < border || @warn preprint * " allocated more memory than usual: " *
		Base.format_bytes(bytes) * " vs. " * Base.format_bytes(border)
end


# Single choice variable
model = :CapitalAdjustModel
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5, :β => 0.9)
p_1_All = createmodel(model; dipar..., intdim = :All)
p_1_Separable = createmodel(model; dipar..., intdim = :Separable)
p_1_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_1_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardcall=>:pre)

sol_1_All = solve(p_1_All; diopt...)
val, t, bytes = @timed solve(p_1_All; diopt...)
warn_time(t, 0.12, "CapitalAdjustModel All: Solver")
warn_memory(bytes, 75555555, "CapitalAdjustModel All: Solver")

sol_1_Separable = solve(p_1_Separable; diopt...)
val, t, bytes = @timed solve(p_1_Separable; diopt...)
warn_time(t, 0.07, "CapitalAdjustModel Separable: Solver")
warn_memory(bytes, 19000000, "CapitalAdjustModel Separable: Solver")

sol_1_Separable_States = solve(p_1_Separable_States; diopt...)
val, t, bytes = @timed solve(p_1_Separable_States; diopt...)
warn_time(t, 0.01, "CapitalAdjustModel Separable_States: Solver")
warn_memory(bytes, 650000, "CapitalAdjustModel Separable_States: Solver")

sol_1_Separable_ExogStates = solve(p_1_Separable_ExogStates; diopt...)
val, t, bytes = @timed solve(p_1_Separable_ExogStates; diopt...)
warn_time(t, 0.01, "CapitalAdjustModel Separable_ExogStates: Solver")
warn_memory(bytes, 210000, "CapitalAdjustModel Separable_ExogStates: Solver")

# SIMULATE
shocks = drawshocks(p_1_All, nFirms=10^4, nPeriods=50)

simulate(sol_1_All, shocks)
val, t, bytes = @timed simulate(sol_1_All, shocks)
warn_time(t, 0.25, "CapitalAdjustModel All: Simulator")
warn_memory(bytes, 221000000, "CapitalAdjustModel All: Simulator")

simulate(sol_1_Separable, shocks)
val, t, bytes = @timed simulate(sol_1_Separable, shocks)
warn_time(t, 0.25, "CapitalAdjustModel Separable: Simulator")
warn_memory(bytes, 260000000, "CapitalAdjustModel Separable: Simulator")

simulate(sol_1_Separable_States, shocks)
val, t, bytes = @timed simulate(sol_1_Separable_States, shocks)
warn_time(t, 0.24, "CapitalAdjustModel Separable_States: Simulator")
warn_memory(bytes, 260000000, "CapitalAdjustModel Separable_States: Simulator")

simulate(sol_1_Separable_ExogStates, shocks)
simulate(sol_1_Separable_ExogStates, shocks)
val, t, bytes = @timed simulate(sol_1_Separable_ExogStates, shocks)
warn_time(t, 0.23, "CapitalAdjustModel Separable_ExogStates: Simulator")
warn_memory(bytes, 190000000, "CapitalAdjustModel Separable_ExogStates: Simulator")


# Two choice variables
model = :CapitalAdjustModel2
dipar = Dict(:nK => 10, :nN => 10, :nz => 5, :β => 0.9)
p_2_All = createmodel(model; dipar..., intdim = :All)
p_2_Sep = createmodel(model; dipar..., intdim = :Separable)
p_2_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_2_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardcall=>:pre)

sol_2_All = solve(p_2_All; diopt...)
val, t, bytes = @timed solve(p_2_All; diopt...)
warn_time(t, 2.1, "CapitalAdjustModel2 All: Solver")
warn_memory(bytes, 860000000, "CapitalAdjustModel2 All: Solver")

sol_2_Sep = solve(p_2_Sep; diopt...)
val, t, bytes = @timed solve(p_2_Sep; diopt...)
warn_time(t, 1., "CapitalAdjustModel2 Separable: Solver")
warn_memory(bytes, 800000000, "CapitalAdjustModel2 Separable: Solver")

sol_2_Sep_States = solve(p_2_Sep_States; diopt...)
val, t, bytes = @timed solve(p_2_Sep_States; diopt...)
warn_time(t, 0.1, "CapitalAdjustModel2 Separable_States: Solver")
warn_memory(bytes, 2000000, "CapitalAdjustModel2 Separable_States: Solver")

sol_2_Sep_ExogStates = solve(p_2_Sep_ExogStates; diopt...)
val, t, bytes = @timed solve(p_2_Sep_ExogStates; diopt...)
warn_time(t, 0.1, "CapitalAdjustModel2 Separable_ExogStates: Solver")
warn_memory(bytes, 552000, "CapitalAdjustModel2 Separable_ExogStates: Solver")

# SIMULATE
shocks = drawshocks(p_2_All, nFirms=10^3, nPeriods=50)

simulate(sol_2_All, shocks)
val, t, bytes = @timed simulate(sol_2_All, shocks)
warn_time(t, 0.12, "CapitalAdjustModel2 All: Simulator")
warn_memory(bytes, 62000000, "CapitalAdjustModel2 All: Simulator")

simulate(sol_2_Sep, shocks)
val, t, bytes = @timed simulate(sol_2_Sep, shocks)
warn_time(t, 0.07, "CapitalAdjustModel2 Separable: Simulator")
warn_memory(bytes, 48000000, "CapitalAdjustModel2 Separable: Simulator")

simulate(sol_2_Sep_States, shocks)
val, t, bytes = @timed simulate(sol_2_Sep_States, shocks)
warn_time(t, 0.07, "CapitalAdjustModel2 Separable_States: Simulator")
warn_memory(bytes, 48000000, "CapitalAdjustModel2 Separable_States: Simulator")

simulate(sol_2_Sep_ExogStates, shocks)
val, t, bytes = @timed simulate(sol_2_Sep_ExogStates, shocks)
warn_time(t, 0.07, "CapitalAdjustModel2 Separable_ExogStates: Simulator")
warn_memory(bytes, 48000000, "CapitalAdjustModel2 Separable_ExogStates: Simulator")
