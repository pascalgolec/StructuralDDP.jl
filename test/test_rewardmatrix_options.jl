p_1 = createmodel(:CapitalAdjustModel; nK = 40, nz = 11, γ = 2., intdim = :Separable)
sol_1_prebuild = solve(p_1; rewardcall=:pre)
sol_1_prebuild_partial = solve(p_1; rewardcall=:pre_partial)
sol_1_nobuild = solve(p_1; rewardcall=:jit)

p_2 = createmodel(:CapitalAdjustModel2; nK = 15, nN = 10, nz = 3, intdim=:Separable_ExogStates,)
diopt = Dict(:concavity => true, :monotonicity => true)
sol_2_prebuild = solve(p_2; diopt..., rewardcall=:pre)
sol_2_prebuild_partial = solve(p_2; diopt..., rewardcall=:pre_partial)
sol_2_nobuild = solve(p_2; diopt..., rewardcall=:jit)

@testset "reward matrix options" begin

	compare(sol_1_nobuild, sol_1_prebuild_partial, str="SingleChoiceVar & jit & pre_partial", tol=mytol)
	compare(sol_1_nobuild, sol_1_prebuild, str="SingleChoiceVar & jit & pre", tol=mytol)
	compare(sol_1_prebuild, sol_1_prebuild_partial, str="SingleChoiceVar & pre & pre_partial", tol=mytol)

	compare(sol_2_nobuild, sol_2_prebuild_partial, str="TwoChoiceVar & jit & pre_partial", tol=mytol)
	compare(sol_2_nobuild, sol_2_prebuild, str="TwoChoiceVar &jit & pre", tol=mytol)
	compare(sol_2_prebuild, sol_2_prebuild_partial, str="TwoChoiceVar &pre & pre_partial", tol=mytol)

end
