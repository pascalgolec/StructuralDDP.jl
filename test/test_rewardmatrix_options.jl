using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4


p_neo = createmodel(:NeoClassicalSimple; nK = 40, nz = 11, γ = 0.5, F = 0., τ = 0.3)
sol_neo_prebuild = solve(p_neo; intdim = :separable, rewardmat=:prebuild)
sol_neo_prebuild_partial = solve(p_neo; intdim = :separable, rewardmat=:prebuild_partial)
sol_neo_nobuild = solve(p_neo; intdim = :separable, rewardmat=:nobuild)

p_intan = createmodel(:Intangible; nK = 25, nN = 20, nz = 3)
optdict = Dict(:intdim=>:separable, :monotonicity=>true, :concavity=>true)
sol_intan_prebuild = solve(p_intan; optdict..., rewardmat=:prebuild)
sol_intan_prebuild_partial = solve(p_intan; optdict..., rewardmat=:prebuild_partial)
sol_intan_nobuild = solve(p_intan; optdict..., rewardmat=:nobuild)

@testset "reward matrix options" begin

	compare_solutions("SingleChoiceVar & nobuild & prebuild_partial", sol_neo_nobuild, sol_neo_prebuild_partial)
	compare_solutions("SingleChoiceVar & nobuild & prebuild", sol_neo_nobuild, sol_neo_prebuild)
	compare_solutions("SingleChoiceVar & prebuild & prebuild_partial", sol_neo_prebuild, sol_neo_prebuild_partial)

	compare_solutions("TwoChoiceVar & nobuild & prebuild_partial", sol_intan_nobuild, sol_intan_prebuild_partial)
	compare_solutions("TwoChoiceVar &nobuild & prebuild", sol_intan_nobuild, sol_intan_prebuild)
	compare_solutions("TwoChoiceVar &prebuild & prebuild_partial", sol_intan_prebuild, sol_intan_prebuild_partial)

end
