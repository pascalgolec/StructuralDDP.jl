struct DDPSimulation{NS,NC,T}
	prob::DDP
	sol::AbstractDDPSolution{NS,NC}
	value::Union{Array{T,3}, Nothing}
	states::Array{T,3}
	policy::Array{T,3}
end
value(sim::DDPSimulation) = sim.value
states(sim::DDPSimulation) = sim.states
policy(sim::DDPSimulation) = sim.policy

"""Construct a DataFrame from the simulation."""
function DataFrame(sim::DDPSimulation)
	(nStateVars, nPeriods, nFirms) = size(states(sim))
	nChoiceVars = size(policy(sim), 1)

	mFirm = repeat(collect(1:nFirms)', nPeriods, 1)
	mPeriod = repeat(collect(0:nPeriods-1), 1, nFirms)

	df = DataFrame()
	df.firm = mFirm[:]
	df.period = mPeriod[:]

	[df[Symbol("state_" * string(i))] = states(sim)[i,:,:][:] for i = 1 : nStateVars]
	[df[Symbol("choice_" * string(i))] = policy(sim)[i,:,:][:] for i = 1 : nChoiceVars]
	if value(sim) != nothing
		df[:value] = value(sim)[1,:,:][:]
	end
	return df
end

"""Construct an `variable x periods x firms` array from the simulation.
The order of the variables is states, choices and value."""
function Array(sim::DDPSimulation)
	if value(sim) == nothing
		return cat(states(sim), policy(sim), dims=1)
	else
		return cat(states(sim), policy(sim), value(sim), dims=1)
	end
end
