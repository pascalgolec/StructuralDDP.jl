
using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 11, :γ => 0.5, :F => 0., :τ => 0.3)
p_neoclassical = createmodel(model; dipar...)

sol_neoclassical = solve(p_neoclassical;
    intdim = :separable, monotonicity = false, concavity = false)
sol_neoclassical_mon = solve(p_neoclassical;
    intdim = :separable, monotonicity = true, concavity = false)
sol_neoclassical_conc = solve(p_neoclassical;
    intdim = :separable, monotonicity = false, concavity = true)

p_intan = createmodel(:Intangible; nK = 25, nN = 20, nz = 3)
optdict = Dict(:intdim=>:separable, :rewardmat=>:prebuild_partial)
sol_intan = solve(p_intan; optdict..., monotonicity=[false,false], concavity=[false,false])
sol_intan_mon = solve(p_intan; optdict..., monotonicity=[true,true], concavity=[false,false])
sol_intan_conc = solve(p_intan; optdict..., monotonicity=[false,false], concavity=[true,true])

@testset "monotonicity concavity" begin

    compare_solutions("Neoclassical_monotonicity", sol_neoclassical, sol_neoclassical_mon)
    compare_solutions("Neoclassical_concavity", sol_neoclassical, sol_neoclassical_conc)
    compare_solutions("Intangible_monotonicity", sol_intan, sol_intan_mon)
    compare_solutions("Intangible_concavity", sol_intan, sol_intan_conc)

end
