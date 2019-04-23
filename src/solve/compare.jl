"""Compare different solutions, useful for when testing monotonicity and concavity."""
function compare(sol1::AbstractDDPSolution,
        sol2::AbstractDDPSolution;
		str::String = "compare solutions",
        tol::Real = 1e-4)
    @testset "$str" begin
        for field in fieldnames(typeof(sol1))
			field1 = getfield(sol1,field)
			field2 = getfield(sol2,field)
            if typeof(field1) <: Tuple
                [@test isapprox(field1[i], field2[i], rtol=tol) for i = 1 : length(field1)]
            elseif typeof(field1) <: Array
                @test isapprox(field1, field2, rtol=tol)
            end
        end
    end
end
