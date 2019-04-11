using DiscreteDynamicProgramming
include("utils.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

# Single choice variable
dipar = Dict(:nK => 40, :nz => 5, :γ => 2., :F => 0., :τ => 0.)
p_neoclassical = createmodel(:NeoClassicalSimple; dipar..., intdim = :Separable_ExogStates)
sol_neoclassical = solve(p_neoclassical; monotonicity = false, concavity = false)
sol_neoclassical_mon = solve(p_neoclassical; monotonicity = true, concavity = false)
sol_neoclassical_conc = solve(p_neoclassical; monotonicity = false, concavity = true)

# two choice variables
dipar = Dict(:nK => 10, :nN => 7, :nz => 3)
p_int = createmodel(:Intangible; dipar..., intdim = :Separable_ExogStates)
optdict = Dict(:rewardmat=>:nobuild)
sol_int = solve(p_int; optdict..., monotonicity=[false,false], concavity=[false,false])
sol_int_mon = solve(p_int; optdict..., monotonicity=[true,true], concavity=[false,false])
sol_int_conc = solve(p_int; optdict..., monotonicity=[false,false], concavity=[true,true])

@testset "monotonicity concavity" begin

    compare_solutions("Neoclassical_monotonicity", sol_neoclassical, sol_neoclassical_mon, mytol)
    compare_solutions("Neoclassical_concavity", sol_neoclassical, sol_neoclassical_conc, mytol)
    compare_solutions("Intangible_monotonicity", sol_int, sol_int_mon, mytol)
    compare_solutions("Intangible_concavity", sol_int, sol_int_conc, mytol)

end
