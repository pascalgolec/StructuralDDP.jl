# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"

inbounds(x::T1, xmin::T2, xmax::T3) where {T1<:Real, T2<:Real, T3<:Real} =
	max(xmin, min(x, xmax))

getiterator(x::NTuple{N,AbstractVector}) where N = Iterators.product(x...)
getiterator(x::NTuple{1,AbstractVector}) = x[1]

const FuncOrNothing = Union{Function,Nothing}

"""Returns a modified version of the function f that if forced to stay between
the bounds defined by the tuple of vectors tVectors"""
function wrapinbounds(f, tVectors::NTuple{N,Vector{C}}) where {N,C}
	lowerbounds::NTuple{N,C} = minimum.(tVectors)
	upperbounds::NTuple{N,C} = maximum.(tVectors)
	return function f2(args...)
		inbounds.(f(args...), lowerbounds, upperbounds)
	end
end
function wrapinbounds(f, tVectors::NTuple{1,Vector{C}}) where {C}
	lowerbounds::C = minimum(tVectors[1])
	upperbounds::C = maximum(tVectors[1])
	return function f1(args...)
		inbounds(f(args...), lowerbounds, upperbounds)
	end
end
