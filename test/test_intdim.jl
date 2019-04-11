
using DiscreteDynamicProgramming
include("utils.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 5, :γ => 0.5, :F => 0., :τ => 0.3)
p_neoclassical_All = createmodel(model; dipar..., intdim = :All)
p_neoclassical_Separable = createmodel(model; dipar..., intdim = :Separable)
p_neoclassical_Separable_States = createmodel(model; dipar..., intdim = :Separable_States)
p_neoclassical_Separable_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

diopt = Dict(:concavity => false, :monotonicity => false, :rewardmat=>:prebuild,
	:initialize_exact=>false)
sol_neoclassical_All = solve(p_neoclassical_All; diopt...)
sol_neoclassical_Separable = solve(p_neoclassical_Separable; diopt...)
sol_neoclassical_Separable_States = solve(p_neoclassical_Separable_States; diopt...)
sol_neoclassical_Separable_ExogStates = solve(p_neoclassical_Separable_ExogStates; diopt...)


model = :Intangible
dipar = Dict(:nK => 10, :nN => 7, :nz => 3)
p_int_All = createmodel(model; dipar..., intdim = :All)
p_int_Sep = createmodel(model; dipar..., intdim = :Separable)
p_int_Sep_States = createmodel(model; dipar..., intdim = :Separable_States)
p_int_Sep_ExogStates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

optdict = Dict(
	:monotonicity=>[false,false], :concavity=>[false,false],
	:rewardmat=>:nobuild, :initialize_exact=>false)
sol_intan_All = solve(p_int_All; optdict...)
sol_intan_Separable = solve(p_int_Sep; optdict...)
sol_intan_Separable_States = solve(p_int_Sep_States; optdict...)
sol_intan_Separable_ExogStates = solve(p_int_Sep_ExogStates; optdict...)


@testset "intdims" begin

	@testset "compare solution methods" begin

        compare_solutions("Neoclassical_1st_2nd", sol_neoclassical_All, sol_neoclassical_Separable, mytol)
        compare_solutions("Neoclassical_2nd_3rd", sol_neoclassical_Separable, sol_neoclassical_Separable_States, mytol)
        compare_solutions("Neoclassical_3rd_4th", sol_neoclassical_Separable_States, sol_neoclassical_Separable_ExogStates, mytol)

		compare_solutions("Intangible_1st_2nd", sol_intan_All, sol_intan_Separable, mytol)
		compare_solutions("Intangible_2nd_3rd", sol_intan_Separable, sol_intan_Separable_States, mytol)
        compare_solutions("Intangible_3rd_4th", sol_intan_Separable_States, sol_intan_Separable_ExogStates, mytol)

    end

    # @testset "SA" begin
    #
    #     @testset "NeoClassical" begin
    #         @testset "policy" begin
    #             @test sol_neoclassical_SA.meshValFun[1] ≈ 11.47959458057959 rtol=mytol
    #             @test sol_neoclassical_SA.meshValFun[end] ≈ 96.85333912427213 rtol=mytol
    #             @test sol_neoclassical_SA.meshPolFun[1] ≈ 2.5896937673964047 rtol=mytol
    #             @test sol_neoclassical_SA.meshPolFun[end] ≈ 51.79387534792809 rtol=mytol
    #         end
    #
    #         @testset "inital policy" begin
    #             @test sol_neoclassical_SA.mV0[1] ≈ 7.336084552745342 rtol=mytol
    #             @test sol_neoclassical_SA.mV0[end] ≈ 22.660418124491528 rtol=mytol
    #             @test sol_neoclassical_SA.mPolFun0[1] ≈ 2.5896937673964047 rtol=mytol
    #             @test sol_neoclassical_SA.mPolFun0[end] ≈ 19.081026821842094 rtol=mytol
    #         end
    #     end
    # end # SA



end # intdims
