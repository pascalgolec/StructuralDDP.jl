inbounds(x::Float64, xmin::Float64, xmax::Float64) = max(xmin, min(x, xmax))

# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"

"""Retreive the state vectors that are not choice vectors from the state vector tuple."""
getnonchoicevars(p::DDM) = getnonchoicevars(p.tStateVectors, p.tChoiceVectors)
function getnonchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC}
	allvars = tuple(1:NS...)
	nonchoicevars = tuple(setdiff(Set(allvars), Set(tchoicevars))...)
	return getindex(tStateVectors, collect(nonchoicevars))
end

getnonchoicevarszero(p::DDM) = getnonchoicevars(p.tStateVectors, p.tChoiceVectorsZero)


"""Retreive the choice vectors from the state vector tuple."""
getchoicevars(p::DDM) = getchoicevars(p.tStateVectors, p.tChoiceVectors)
getchoicevars(tStateVectors::NTuple{NS,T}, tchoicevars::NTuple{NC, Int64}) where {NS, T, NC} =
	getindex(tStateVectors, collect(tchoicevars))
"""Retreive nothing if the choice vectors are provided."""
getchoicevars(tStateVectors::NTuple{N1,T1}, tChoiceVectors::NTuple{N2,Vector{T2}}) where
	{N1, T1, N2, T2} = tChoiceVectors

"""Retreive the choice vectors to find policy at t=0."""
getchoicevarszero(p::DDM) = getchoicevars(p.tStateVectors, p.tChoiceVectorsZero)
