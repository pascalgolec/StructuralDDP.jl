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
