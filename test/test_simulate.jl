using DataFrames

dipar = Dict(:nK => 40, :nz => 5, :γ => 2.)
p_1_all = createmodel(:CapitalAdjustModel; dipar..., :intdim=>:All)
p_1_sep = createmodel(:CapitalAdjustModel; dipar..., :intdim=>:Separable)
p_1_sep_states = createmodel(:CapitalAdjustModel; dipar..., :intdim=>:Separable_States)
p_1_sep_exogstates = createmodel(:CapitalAdjustModel; dipar..., :intdim=>:Separable_ExogStates)

optdict = Dict(:monotonicity=>true, :concavity=>true,
	:rewardcall=>:pre_partial)


sol_1_all = solve(p_1_all; optdict...)
sol_1_sep = solve(p_1_sep; optdict...)
sol_1_sep_states = solve(p_1_sep_states; optdict...)
sol_1_sep_exogstates = solve(p_1_sep_exogstates; optdict...)

shocks_1 = drawshocks(p_1_all, nPeriods=40, nFirms=100)

@testset "Simulate CapitalAdjustModel" begin
	@testset "conversion" begin
		sim = simulate(sol_1_sep_exogstates, shocks_1, get_value=true)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
		sim = simulate(sol_1_sep_exogstates, shocks_1, get_value=false)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
	end
	@testset "initialize_inexact" begin
		sim_opt = Dict(:initialize_exact=>false, :get_value=>true)
	    sim_1 = simulate(sol_1_all, shocks_1; sim_opt...)
	    sim_2 = simulate(sol_1_sep, shocks_1; sim_opt...)
	    sim_3 = simulate(sol_1_sep_states, shocks_1; sim_opt...)
	    sim_4 = simulate(sol_1_sep_exogstates, shocks_1; sim_opt...)
		@test isapprox(Array(sim_1), Array(sim_2), rtol=mytol)
		@test isapprox(Array(sim_2), Array(sim_3), rtol=mytol)
		@test isapprox(Array(sim_3), Array(sim_4), rtol=mytol)
	end
	@testset "initialize_exact" begin
		sim_opt = Dict(:initialize_exact=>true, :get_value=>true)
	    sim_1 = simulate(sol_1_all, shocks_1; sim_opt...)
	    sim_2 = simulate(sol_1_sep, shocks_1; sim_opt...)
	    sim_3 = simulate(sol_1_sep_states, shocks_1; sim_opt...)
	    sim_4 = simulate(sol_1_sep_exogstates, shocks_1; sim_opt...)
		@test isapprox(Array(sim_1), Array(sim_2), rtol=mytol)
		@test isapprox(Array(sim_2), Array(sim_3), rtol=mytol)
		@test isapprox(Array(sim_3), Array(sim_4), rtol=mytol)
	end
end


model = :CapitalAdjustModel2
dipar = Dict(:nK => 15, :nN => 10, :nz => 5, :ρ => 0.6)
p_2_all = createmodel(model; dipar..., :intdim=>:All)
p_2_sep = createmodel(model; dipar..., :intdim=>:Separable)
p_2_sep_states = createmodel(model; dipar..., :intdim=>:Separable_States)
p_2_sep_exogstates = createmodel(model; dipar..., :intdim=>:Separable_ExogStates)

optdict = Dict(:monotonicity=>[true,true], :concavity=>[true,true],
	:rewardcall=>:pre_partial)
optdict = Dict(:monotonicity=>false, :concavity=>false,
	:rewardcall=>:pre_partial)

sol_2_all = solve(p_2_all; optdict...)
sol_2_sep = solve(p_2_sep; optdict...)
sol_2_sep_states = solve(p_2_sep_states; optdict...)
sol_2_sep_exogstates = solve(p_2_sep_exogstates; optdict...)

shocks_2 = drawshocks(p_2_all, nPeriods=20, nFirms=10)

isapprox(sol_2_all, sol_2_sep)
isapprox(sol_2_all, sol_2_sep_states)
isapprox(sol_2_all, sol_2_sep_exogstates)
@testset "Simulate CapitalAdjustModel2" begin
	@testset "conversion" begin
		sim = simulate(sol_2_sep_exogstates, shocks_2, get_value=true)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
		sim = simulate(sol_2_sep_exogstates, shocks_2, get_value=false)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
	end
	@testset "initialize_inexact" begin
		sim_opt = Dict(:initialize_exact=>false, :get_value=>true)
	    sim_1 = simulate(sol_2_all, shocks_2; sim_opt...)
	    sim_2 = simulate(sol_2_sep, shocks_2; sim_opt...)
	    sim_3 = simulate(sol_2_sep_states, shocks_2; sim_opt...)
	    sim_4 = simulate(sol_2_sep_exogstates, shocks_2; sim_opt...)
		@test isapprox(Array(sim_1), Array(sim_2), rtol=mytol)
		@test isapprox(Array(sim_2), Array(sim_3), rtol=mytol)
		@test isapprox(Array(sim_3), Array(sim_4), rtol=mytol)
	end
	@testset "initialize_exact" begin
		sim_opt = Dict(:initialize_exact=>true, :get_value=>true)
	    sim_1 = simulate(sol_2_all, shocks_2; sim_opt...)
	    sim_2 = simulate(sol_2_sep, shocks_2; sim_opt...)
	    sim_3 = simulate(sol_2_sep_states, shocks_2; sim_opt...)
	    sim_4 = simulate(sol_2_sep_exogstates, shocks_2; sim_opt...)
		@test isapprox(Array(sim_1), Array(sim_2), rtol=mytol)
		@test isapprox(Array(sim_2), Array(sim_3), rtol=mytol)
		@test isapprox(Array(sim_3), Array(sim_4), rtol=mytol)
	end
end
