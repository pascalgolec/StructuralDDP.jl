
using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 11, :γ => 0.5, :F => 0., :τ => 0.3)
p_neoclassical = createmodel(model; dipar...)

diopt = Dict(:concavity => false, :monotonicity => false,
	:rewardmat=>:prebuild_partial)
sol_neoclassical_SA = solve(p_neoclassical, intdim = :SA)
sol_neoclassical_separable = solve(p_neoclassical; diopt..., intdim = :separable)
sol_neoclassical_intermediate = solve(p_neoclassical; diopt..., intdim = :intermediate)

ptest = createmodel(:Intangible; nK = 15, nN = 10, nz = 3)
optdict = Dict(
	:monotonicity=>[true,true],
	:concavity=>[true,true],
	:rewardmat=>:prebuild_partial,
	)
sol_intan_SA = solve(ptest; intdim = :SA)
sol_intan_separable = solve(ptest; optdict..., intdim = :separable)
sol_intan_intermediate = solve(ptest; optdict..., intdim = :intermediate)

@testset "intdims" begin

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

    @testset "compare solution methods" begin

        compare_solutions("Neoclassical_SA_sep", sol_neoclassical_SA, sol_neoclassical_separable)
        compare_solutions("Neoclassical_SA_int", sol_neoclassical_SA, sol_neoclassical_intermediate)
        compare_solutions("Neoclassical_sep_int", sol_neoclassical_separable, sol_neoclassical_intermediate)

		compare_solutions("Intangible_SA_sep", sol_intan_SA, sol_intan_separable)
        compare_solutions("Intangible_SA_int", sol_intan_SA, sol_intan_intermediate)
        compare_solutions("Intangible_sep_int", sol_intan_separable, sol_intan_intermediate)

    end

end # intdims
