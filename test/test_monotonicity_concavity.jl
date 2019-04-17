using DiscreteDynamicProgramming
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

    compare(sol_neoclassical, sol_neoclassical_mon, str="Neoclassical_monotonicity",  tol=mytol)
    compare(sol_neoclassical, sol_neoclassical_conc, str="Neoclassical_concavity", tol=mytol)
    compare(sol_int, sol_int_mon, str="Intangible_monotonicity", tol=mytol)
    compare(sol_int, sol_int_conc, str="Intangible_concavity", tol=mytol)

end
