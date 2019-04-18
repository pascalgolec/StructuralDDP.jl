using DiscreteDynamicProgramming

function warn_time(t, border, preprint)
	t < border || @warn preprint * "Solver took longer than usual: $t vs $border seconds"
end
function warn_memory(bytes, border, preprint)
	bytes < border || @warn preprint * "Solver allocated more memory than usual: " *
		Base.format_bytes(bytes) * " vs. " * Base.format_bytes(border)
end


# Single choice variable
model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5, :F => 0., :τ => 0.3, :β => 0.9)
p_neoclassical_All = createmodel(model; dipar..., intdim = :All)
p_neoclassical_Separable = createmodel(model; dipar..., intdim = :Separable)
p_neoclassical_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_neoclassical_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardcall=>:pre,
	:initialize_exact=>false)

solve(p_neoclassical_All; diopt...)
val, t, bytes = @timed solve(p_neoclassical_All; diopt...)
warn_time(t, 0.12, "Neo All: ")
warn_memory(bytes, 75555555, "Neo All: ")

solve(p_neoclassical_Separable; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable; diopt...)
warn_time(t, 0.07, "Neo Separable: ")
warn_memory(bytes, 19000000, "Neo Separable: ")

sol_neoclassical_Separable_States = solve(p_neoclassical_Separable_States; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable_States; diopt...)
warn_time(t, 0.01, "Neo Separable_States: ")
warn_memory(bytes, 650000, "Neo Separable_States: ")

sol_neoclassical_Separable_ExogStates = solve(p_neoclassical_Separable_ExogStates; diopt...)
val, t, bytes = @timed solve(p_neoclassical_Separable_ExogStates; diopt...)
warn_time(t, 0.01, "Neo Separable_ExogStates: ")
warn_memory(bytes, 200000, "Neo Separable_ExogStates: ")


# Two choice variables
model = :Intangible
dipar = Dict(:nK => 10, :nN => 10, :nz => 5, :β => 0.9)
p_int_All = createmodel(model; dipar..., intdim = :All)
p_int_Sep = createmodel(model; dipar..., intdim = :Separable)
p_int_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_int_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

optdict = Dict(
	:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardcall=>:jit, :initialize_exact=>false)

val, t, bytes = @timed solve(p_int_All; diopt...)
val, t, bytes = @timed solve(p_int_All; diopt...)
val, t, bytes = @timed solve(p_int_All; diopt...)
warn_time(t, 2., "Intangible All: ")
warn_memory(bytes, 800000000, "Intangible All: ")

val, t, bytes = @timed solve(p_int_Sep; diopt...)
val, t, bytes = @timed solve(p_int_Sep; diopt...)
val, t, bytes = @timed solve(p_int_Sep; diopt...)
warn_time(t, 1., "Intangible Separable: ")
warn_memory(bytes, 800000000, "Intangible Separable: ")

val, t, bytes = @timed solve(p_int_Sep_States; diopt...)
val, t, bytes = @timed solve(p_int_Sep_States; diopt...)
val, t, bytes = @timed solve(p_int_Sep_States; diopt...)
warn_time(t, 0.1, "Intangible Separable_States: ")
warn_memory(bytes, 2000000, "Intangible Separable_States: ")

val, t, bytes = @timed solve(p_int_Sep_ExogStates; diopt...)
val, t, bytes = @timed solve(p_int_Sep_ExogStates; diopt...)
val, t, bytes = @timed solve(p_int_Sep_ExogStates; diopt...)
warn_time(t, 0.1, "Intangible Separable_ExogStates: ")
warn_memory(bytes, 550000, "Intangible Separable_ExogStates: ")
