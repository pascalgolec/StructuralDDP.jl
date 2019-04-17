inbounds(x::T1, xmin::T2, xmax::T3) where {T1<:Real, T2<:Real, T3<:Real} =
	max(xmin, min(x, xmax))

# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"

"""Compare different solutions, useful for when testing monotonicity and concavity."""
function compare(sol1::DDPSolution,
        sol2::DDPSolution;
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
