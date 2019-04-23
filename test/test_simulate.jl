using DiscreteDynamicProgramming
using Test
using DataFrames

# todo: find a better way to test rather than just calling the simulate function

dipar = Dict(:nK => 40, :nz => 5, :γ => 2., :F => 0., :τ => 0.3)
p_neo_all = createmodel(:NeoClassicalSimple; dipar..., :intdim=>:All)
p_neo_sep = createmodel(:NeoClassicalSimple; dipar..., :intdim=>:Separable)
p_neo_sep_states = createmodel(:NeoClassicalSimple; dipar..., :intdim=>:Separable_States)
p_neo_sep_exogstates = createmodel(:NeoClassicalSimple; dipar..., :intdim=>:Separable_ExogStates)

optdict = Dict(
	:monotonicity=>true,
	:concavity=>true,
	:rewardcall=>:pre_partial)


sol_neo_all = solve(p_neo_all; optdict...)
sol_neo_sep = solve(p_neo_sep; optdict...)
sol_neo_sep_states = solve(p_neo_sep_states; optdict...)
sol_neo_sep_exogstates = solve(p_neo_sep_exogstates; optdict...)

shocks_neo = drawshocks(p_neo_all, nPeriods=40, nFirms=100)

@testset "Simulate Neoclassical" begin
	@testset "initialize_inexact" begin
	    simulate(sol_neo_all, shocks_neo, initialize_exact=false)
	    simulate(sol_neo_sep, shocks_neo, initialize_exact=false)
	    simulate(sol_neo_sep_states, shocks_neo, initialize_exact=false)
	    simulate(sol_neo_sep_exogstates, shocks_neo, initialize_exact=false)
	end
	@testset "initialize_exact" begin
	    simulate(sol_neo_all, shocks_neo, initialize_exact=true)
	    simulate(sol_neo_sep, shocks_neo, initialize_exact=true)
	    simulate(sol_neo_sep_states, shocks_neo, initialize_exact=true)
	    simulate(sol_neo_sep_exogstates, shocks_neo, initialize_exact=true)
	end
	@testset "conversion" begin
		sim = simulate(sol_neo_sep_exogstates, shocks_neo, get_value=true)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
		sim = simulate(sol_neo_sep_exogstates, shocks_neo, get_value=false)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
	end
end


model = :Intangible
dipar = Dict(:nK => 15, :nN => 10, :nz => 5, :τ => 0.0, :ρ => 0.6)
p_int_all = createmodel(model; dipar..., :intdim=>:All)
p_int_sep = createmodel(model; dipar..., :intdim=>:Separable)
p_int_sep_states = createmodel(model; dipar..., :intdim=>:Separable_States)
p_int_sep_exogstates = createmodel(model; dipar..., :intdim=>:Separable_ExogStates)

optdict = Dict(
	:monotonicity=>true,
	:concavity=>true,
	:rewardcall=>:pre_partial)

sol_int_all = solve(p_int_all; optdict...)
sol_int_sep = solve(p_int_sep; optdict...)
sol_int_sep_states = solve(p_int_sep_states; optdict...)
sol_int_sep_exogstates = solve(p_int_sep_exogstates; optdict...)

shocks_int = drawshocks(p_int_all, nPeriods=20, nFirms=10)

@testset "Simulate Intangible" begin
	@testset "initialize_inexact" begin
	    simulate(sol_int_all, shocks_int, initialize_exact=false)
	    simulate(sol_int_sep, shocks_int, initialize_exact=false)
	    simulate(sol_int_sep_states, shocks_int, initialize_exact=false)
	    simulate(sol_int_sep_exogstates, shocks_int, initialize_exact=false)
	end
	@testset "initialize_exact" begin
	    simulate(sol_int_all, shocks_int, initialize_exact=true)
	    simulate(sol_int_sep, shocks_int, initialize_exact=true)
	    simulate(sol_int_sep_states, shocks_int, initialize_exact=true)
	    simulate(sol_int_sep_exogstates, shocks_int, initialize_exact=true)
	end
	@testset "conversion" begin
		sim = simulate(sol_int_sep_exogstates, shocks_neo, get_value=true)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
		sim = simulate(sol_int_sep_exogstates, shocks_neo, get_value=false)
	    @test typeof(DataFrame(sim)) <: DataFrame
		@test typeof(Array(sim)) <: AbstractArray
	end
end
