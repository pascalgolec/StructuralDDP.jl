# Single choice variable
dipar = Dict(:nK => 40, :nz => 5, :γ => 2., :F=>0.02)
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

    @test isapprox(sol_1, sol_1_mon,  rtol=mytol)
    @test isapprox(sol_1, sol_1_conc, rtol=mytol)
    @test isapprox(sol_2, sol_2_mon, rtol=mytol)
    @test isapprox(sol_2, sol_2_conc, rtol=mytol)

end
