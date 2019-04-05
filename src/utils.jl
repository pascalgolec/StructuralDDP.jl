inbounds(x::Float64, xmin::Float64, xmax::Float64) = max(xmin, min(x, xmax))

# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"

function getnonchoicevars(tStateVectors::NTuple{n,T}, tchoicevars::Tuple{Int64}) where {n, T}
	allvars = tuple(1:n...)
	nonchoicevars = tuple(setdiff(Set(allvars), Set(tchoicevars))...)
	return getindex(tStateVectors, collect(nonchoicevars))
end

getchoicevars(tStateVectors::NTuple{n,T}, tchoicevars::Tuple{Int64}) where {n, T} =
	getindex(tStateVectors, collect(tchoicevars))
