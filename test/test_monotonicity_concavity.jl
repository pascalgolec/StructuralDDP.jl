# Single choice variable
dipar = Dict(:nK => 40, :nz => 5, :γ => 2.)
p_1 = createmodel(:CapitalAdjustModel; dipar..., intdim = :Separable_ExogStates)
sol_1 = solve(p_1; monotonicity = false, concavity = false)
sol_1_mon = solve(p_1; monotonicity = true, concavity = false)
sol_1_conc = solve(p_1; monotonicity = false, concavity = true)

# two choice variables
dipar = Dict(:nK => 10, :nN => 7, :nz => 3)
p_2 = createmodel(:CapitalAdjustModel2; dipar..., intdim = :Separable_ExogStates)
optdict = Dict(:rewardcall=>:jit)
sol_2 = solve(p_2; optdict..., monotonicity=[false,false], concavity=[false,false])
sol_2_mon = solve(p_2; optdict..., monotonicity=[true,true], concavity=[false,false])
sol_2_conc = solve(p_2; optdict..., monotonicity=[false,false], concavity=[true,true])

@testset "monotonicity concavity" begin

    compare(sol_1, sol_1_mon, str="CapitalAdjustModel_monotonicity",  tol=mytol)
    compare(sol_1, sol_1_conc, str="CapitalAdjustModel_concavity", tol=mytol)
    compare(sol_2, sol_2_mon, str="CapitalAdjustModel2_monotonicity", tol=mytol)
    compare(sol_2, sol_2_conc, str="CapitalAdjustModel2_concavity", tol=mytol)

end
