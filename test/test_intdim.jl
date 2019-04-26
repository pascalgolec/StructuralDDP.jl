model = :CapitalAdjustModel
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5)
p_1_All = createmodel(model; dipar..., intdim = :All)
p_1_Separable = createmodel(model; dipar..., intdim = :Separable)
p_1_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_1_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardcall=>:pre)
sol_1_All = solve(p_1_All; diopt...)
sol_1_Separable = solve(p_1_Separable; diopt...)
sol_1_Separable_States = solve(p_1_Separable_States; diopt...)
sol_1_Separable_ExogStates = solve(p_1_Separable_ExogStates; diopt...)


model = :CapitalAdjustModel2
dipar = Dict(:nK => 10, :nN => 7, :nz => 3)
p_2_All = createmodel(model; dipar..., intdim = :All)
p_2_Sep = createmodel(model; dipar..., intdim = :Separable)
p_2_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_2_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

optdict = Dict(
	:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardcall=>:jit)
sol_2_All = solve(p_2_All; optdict...)
sol_2_Separable = solve(p_2_Sep; optdict...)
sol_2_Separable_States = solve(p_2_Sep_States; optdict...)
sol_2_Separable_ExogStates = solve(p_2_Sep_ExogStates; optdict...)


@testset "intdims" begin

	@testset "compare solution methods" begin

        @test isapprox(sol_1_All, sol_1_Separable, rtol=mytol)
        @test isapprox(sol_1_Separable, sol_1_Separable_States, rtol=mytol)
        @test isapprox(sol_1_Separable_States, sol_1_Separable_ExogStates, rtol=mytol)

		@test isapprox(sol_2_All, sol_2_Separable, rtol=mytol)
		@test isapprox(sol_2_Separable, sol_2_Separable_States, rtol=mytol)
        @test isapprox(sol_2_Separable_States, sol_2_Separable_ExogStates,rtol=mytol)

    end

end # intdims
