using Test

function compare_solutions(modelname::String,
        sol1::DiscreteDynamicProgramming.DDPSolution,
        sol2::DiscreteDynamicProgramming.DDPSolution,
        mytol)
    @testset "$modelname" begin
        for field in fieldnames(typeof(sol1))
            # @test getfield(sol1,field)[1] ≈ getfield(sol2,field)[1] rtol=mytol
            # @test getfield(sol1,field)[end] ≈ getfield(sol2,field)[end] rtol=mytol
            if typeof(getfield(sol1,field)) <: Tuple
                @test isapprox(getfield(sol1,field)[1], getfield(sol2,field)[1], rtol=mytol)
            else
                @test isapprox(getfield(sol1,field), getfield(sol2,field), rtol=mytol)
            end
        end
    end
end
