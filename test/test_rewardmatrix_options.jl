using DiscreteDynamicProgramming
using Test
mytol = 1e-4


p_neo = createmodel(:NeoClassicalSimple; nK = 40, nz = 11, γ = 2., F = 0., τ = 0., κ=0.,
	intdim = :Separable)
sol_neo_prebuild = solve(p_neo; rewardmat=:prebuild)
sol_neo_prebuild_partial = solve(p_neo; rewardmat=:prebuild_partial)
sol_neo_nobuild = solve(p_neo; rewardmat=:nobuild)

p_intan = createmodel(:Intangible; nK = 15, nN = 10, nz = 3, intdim=:Separable_ExogStates,)
diopt = Dict(:concavity => true, :monotonicity => true, :initialize_exact=>false)
sol_intan_prebuild = solve(p_intan; diopt..., rewardmat=:prebuild)
sol_intan_prebuild_partial = solve(p_intan; diopt..., rewardmat=:prebuild_partial)
sol_intan_nobuild = solve(p_intan; diopt..., rewardmat=:nobuild)

@testset "reward matrix options" begin

	compare(sol_neo_nobuild, sol_neo_prebuild_partial, str="SingleChoiceVar & nobuild & prebuild_partial", tol=mytol)
	compare(sol_neo_nobuild, sol_neo_prebuild, str="SingleChoiceVar & nobuild & prebuild", tol=mytol)
	compare(sol_neo_prebuild, sol_neo_prebuild_partial, str="SingleChoiceVar & prebuild & prebuild_partial", tol=mytol)

	compare(sol_intan_nobuild, sol_intan_prebuild_partial, str="TwoChoiceVar & nobuild & prebuild_partial", tol=mytol)
	compare(sol_intan_nobuild, sol_intan_prebuild, str="TwoChoiceVar &nobuild & prebuild", tol=mytol)
	compare(sol_intan_prebuild, sol_intan_prebuild_partial, str="TwoChoiceVar &prebuild & prebuild_partial", tol=mytol)

end
