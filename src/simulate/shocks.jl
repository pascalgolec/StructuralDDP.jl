# Note: shocks should be independent of parameters, because in SMM will loop over
# parameters, keeping the shocks fixed

abstract type AbstractDDPShocks end

"""Type that contains the shock draw."""
struct DDPShocks <: AbstractDDPShocks

	"""Dimension of aSim is (dimShocks, nPeriods, nFirms)"""
	aSim::Array{Float64,3}
end

"""If `initialize_exact = true`, then also need shocks at t=0."""
struct DDPShocksZero <: AbstractDDPShocks

	"""Dimension of aSim is (dimShocks, nPeriods, nFirms)"""
	aSim::Array{Float64,3}

	"""Dimension of aSim is (dimShocks, 1, nFirms)"""
	aInit::Array{Float64,3}
end

function drawshocks(p::DDP; nPeriods::Int64 = p.params.nPeriods,
		nFirms::Int64 = p.params.nFirms)

	dimShocks = length(p.shockdist)

	aSim = zeros((dimShocks, nPeriods, nFirms))
	for t=1:nPeriods, f=1:nFirms
		aSim[:,t,f] .= rand(p.shockdist)
	end

	if p.options.initialize == nothing
		return DDPShocks(aSim)
	else
		aInit = zeros((dimShocks, 1, nFirms))
		for f = 1:nFirms
			aInit[:,1,f] .= rand(p.shockdist)
		end
		return DDPShocksZero(aSim, aInit)
	end
end
