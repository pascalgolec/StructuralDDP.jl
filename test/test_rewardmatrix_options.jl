using DiscreteDynamicProgramming
include("utils.jl")
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

	compare_solutions("SingleChoiceVar & nobuild & prebuild_partial", sol_neo_nobuild, sol_neo_prebuild_partial, mytol)
	compare_solutions("SingleChoiceVar & nobuild & prebuild", sol_neo_nobuild, sol_neo_prebuild, mytol)
	compare_solutions("SingleChoiceVar & prebuild & prebuild_partial", sol_neo_prebuild, sol_neo_prebuild_partial, mytol)

	compare_solutions("TwoChoiceVar & nobuild & prebuild_partial", sol_intan_nobuild, sol_intan_prebuild_partial, mytol)
	compare_solutions("TwoChoiceVar &nobuild & prebuild", sol_intan_nobuild, sol_intan_prebuild, mytol)
	compare_solutions("TwoChoiceVar &prebuild & prebuild_partial", sol_intan_prebuild, sol_intan_prebuild_partial, mytol)

end
