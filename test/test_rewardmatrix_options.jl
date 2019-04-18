using DiscreteDynamicProgramming
using Test
mytol = 1e-4


p_neo = createmodel(:NeoClassicalSimple; nK = 40, nz = 11, γ = 2., F = 0., τ = 0., κ=0.,
	intdim = :Separable)
sol_neo_prebuild = solve(p_neo; rewardcall=:pre)
sol_neo_prebuild_partial = solve(p_neo; rewardcall=:pre_partial)
sol_neo_nobuild = solve(p_neo; rewardcall=:jit)

p_intan = createmodel(:Intangible; nK = 15, nN = 10, nz = 3, intdim=:Separable_ExogStates,)
diopt = Dict(:concavity => true, :monotonicity => true, :initialize_exact=>false)
sol_intan_prebuild = solve(p_intan; diopt..., rewardcall=:pre)
sol_intan_prebuild_partial = solve(p_intan; diopt..., rewardcall=:pre_partial)
sol_intan_nobuild = solve(p_intan; diopt..., rewardcall=:jit)

@testset "reward matrix options" begin

	compare(sol_neo_nobuild, sol_neo_prebuild_partial, str="SingleChoiceVar & jit & pre_partial", tol=mytol)
	compare(sol_neo_nobuild, sol_neo_prebuild, str="SingleChoiceVar & jit & pre", tol=mytol)
	compare(sol_neo_prebuild, sol_neo_prebuild_partial, str="SingleChoiceVar & pre & pre_partial", tol=mytol)

	compare(sol_intan_nobuild, sol_intan_prebuild_partial, str="TwoChoiceVar & jit & pre_partial", tol=mytol)
	compare(sol_intan_nobuild, sol_intan_prebuild, str="TwoChoiceVar &jit & pre", tol=mytol)
	compare(sol_intan_prebuild, sol_intan_prebuild_partial, str="TwoChoiceVar &pre & pre_partial", tol=mytol)

end
