inbounds(x::T1, xmin::T2, xmax::T3) where {T1<:Real, T2<:Real, T3<:Real} =
	max(xmin, min(x, xmax))

# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"
