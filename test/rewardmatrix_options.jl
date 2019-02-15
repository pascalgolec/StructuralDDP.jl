# include("../src/DiscreteDynamicProgramming.jl")
using DiscreteDynamicProgramming
# DiscreteDynamicProgramming.DDPSolution

include("helpfunctions.jl")
using Test
mytol = 1e-4
disp = false

ptest = createmodel(NeoClassicalSimple; intdim=:separable, rewardmat=:prebuild,
	nK=100, nz=5, F=0.)
sol1 = solve(ptest, disp=false)

ptest = createmodel(NeoClassicalSimple; intdim=:separable, rewardmat=:prebuild_partial,
	nK=100, nz=5, F=0.)
sol2 = solve(ptest, disp=false)

ptest = createmodel(NeoClassicalSimple; intdim=:separable, rewardmat=:nobuild,
	nK=100, nz=5, F=0.)
sol3 = solve(ptest, disp=false)

compare_solutions("prebuild & partial", sol1, sol2)
compare_solutions("prebuild & nobuild", sol1, sol3)
compare_solutions("prebuild_partial & nobuild", sol2, sol3)
