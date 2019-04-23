"""The integration dimension determines the in- and outputs of the transition function."""
abstract type IntegrationDimension end
const IntDim = IntegrationDimension
abstract type All <: IntDim end
abstract type Separable <: IntDim end
abstract type Separable_States <: IntDim end
abstract type Separable_ExogStates <: IntDim end
const Separable_Union = Union{Separable, Separable_States, Separable_ExogStates}

"""Transition function depending on parameters for dispatch."""
struct Transition{ID<:IntDim,nChoices,nShocks,nExogStates,F}
	f::F
end
(p::Transition)(args...) = p.f(args...)

function get_transition(transfunc::Transition, vStates, args...)
	vStatesPrime = similar(vStates)
	get_transition!(transfunc, vStatesPrime, vStates, args...)
	return vStatesPrime
end

function get_transition!(transfunc::Transition{All}, vStatesPrime,
	vStates, vChoices, vShocks)
	vStatesPrime .= transfunc(vStates, vChoices, vShocks)
end

function get_transition!(transfunc::Transition{Separable,1}, vStatesPrime,
		vStates, vChoices, vShocks)
	vStatesPrime[1] = vChoices
	vStatesPrime[2:end] .= transfunc(vStates, vChoices, vShocks)
end
function get_transition!(transfunc::Transition{Separable,NC}, vStatesPrime,
		vStates, vChoices, vShocks) where NC
	vStatesPrime[1:NC] .= vChoices
	vStatesPrime[1+NC:end] .= transfunc(vStates, vChoices, vShocks)
end
function get_transition!(transfunc::Transition{Separable_States,1}, vStatesPrime,
		vStates, vChoices, vShocks)
	vStatesPrime[1] = vChoices
	vStatesPrime[2:end] .= transfunc(vStates, vShocks)
end
function get_transition!(transfunc::Transition{Separable_States,NC}, vStatesPrime,
		vStates, vChoices, vShocks) where NC
	vStatesPrime[1:NC] .= vChoices
	vStatesPrime[1+NC:end] .= transfunc(vStates, vShocks)
end

function get_transition!(transfunc::Transition{Separable_ExogStates,1,NSh,NSE}, vStatesPrime,
		vStates, vChoices, vShocks) where {NSh,NSE}
	vStatesPrime[1] = vChoices
	if NSE == 1
		vStatesPrime[2] = transfunc(vStates[2], vShocks)
	else
		vStatesPrime[2:end] .= transfunc(vStates[2:end], vShocks)
	end
end
function get_transition!(transfunc::Transition{Separable_ExogStates,NC,NSh,NSE}, vStatesPrime,
		vStates, vChoices, vShocks) where {NC,NSh,NSE}
	vStatesPrime[1:NC] .= vChoices
	if NSE == 1
		vStatesPrime[1+NC:end] .= transfunc(vStates[1+NC], vShocks)
	else
		vSim_it1[1+NC:end] .= transfunc(vSim_it[1+NC:end], mShocks_it)
	end
end
