
using DiscreteDynamicProgramming
# include("utils.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5, :F => 0., :τ => 0.3)
p_neoclassical_All = createmodel(model; dipar..., intdim = :All)
p_neoclassical_Separable = createmodel(model; dipar..., intdim = :Separable)
p_neoclassical_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_neoclassical_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardcall=>:pre,
	:initialize_exact=>false)
sol_neoclassical_All = solve(p_neoclassical_All; diopt...)
sol_neoclassical_Separable = solve(p_neoclassical_Separable; diopt...)
sol_neoclassical_Separable_States = solve(p_neoclassical_Separable_States; diopt...)
sol_neoclassical_Separable_ExogStates = solve(p_neoclassical_Separable_ExogStates; diopt...)


model = :Intangible
dipar = Dict(:nK => 10, :nN => 7, :nz => 3)
p_int_All = createmodel(model; dipar..., intdim = :All)
p_int_Sep = createmodel(model; dipar..., intdim = :Separable)
p_int_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_int_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

optdict = Dict(
	:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardcall=>:jit, :initialize_exact=>false)
sol_intan_All = solve(p_int_All; optdict...)
sol_intan_Separable = solve(p_int_Sep; optdict...)
sol_intan_Separable_States = solve(p_int_Sep_States; optdict...)
sol_intan_Separable_ExogStates = solve(p_int_Sep_ExogStates; optdict...)


@testset "intdims" begin

	@testset "compare solution methods" begin

        compare(sol_neoclassical_All, sol_neoclassical_Separable, str="Neoclassical_1st_2nd", tol=mytol)
        compare(sol_neoclassical_Separable, sol_neoclassical_Separable_States, str="Neoclassical_2nd_3rd", tol=mytol)
        compare(sol_neoclassical_Separable_States, sol_neoclassical_Separable_ExogStates, str="Neoclassical_3rd_4th", tol=mytol)

		compare(sol_intan_All, sol_intan_Separable, str="Intangible_1st_2nd", tol=mytol)
		compare(sol_intan_Separable, sol_intan_Separable_States, str="Intangible_2nd_3rd", tol=mytol)
        compare(sol_intan_Separable_States, sol_intan_Separable_ExogStates, str="Intangible_3rd_4th", tol=mytol)

    end

end # intdims
