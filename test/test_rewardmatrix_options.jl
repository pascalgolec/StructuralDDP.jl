using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4


dipar = Dict(:nK => 40, :nz => 11, :γ => 0.5, :F => 0., :τ => 0.3,
    :intdim => :separable)

p_prebuild = createmodel(:NeoClassicalSimple; dipar..., rewardmat=:prebuild)
sol_prebuild = solve(p_prebuild)

p_prebuild_partial = createmodel(:NeoClassicalSimple;  dipar..., rewardmat=:prebuild_partial)
sol_prebuild_partial = solve(p_prebuild_partial)

p_nobuild = createmodel(:NeoClassicalSimple; dipar..., rewardmat=:nobuild)
sol_nobuild = solve(p_nobuild)

@testset "reward matrix options" begin

	compare_solutions("nobuild & prebuild_partial", sol_nobuild, sol_prebuild_partial)
	compare_solutions("nobuild & prebuild", sol_nobuild, sol_prebuild)
	compare_solutions("prebuild & prebuild_partial", sol_prebuild, sol_prebuild_partial)

end
