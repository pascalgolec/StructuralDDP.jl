
using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 11, :γ => 0.5, :F => 0., :τ => 0.3,
    :concavity => false, :monotonicity => false)
p_neoclassical_SA = createmodel(model; dipar..., intdim = :SA)
p_neoclassical_separable = createmodel(model; dipar..., intdim = :separable)
p_neoclassical_intermediate = createmodel(model; dipar..., intdim = :intermediate)

sol_neoclassical_SA = solve(p_neoclassical_SA)
sol_neoclassical_separable = solve(p_neoclassical_separable)
sol_neoclassical_intermediate = solve(p_neoclassical_intermediate)

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

    end

end # intdims
