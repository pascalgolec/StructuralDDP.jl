# Note: shocks should be independent of parameters, because in SMM will loop over
# parameters, keeping the shocks fixed

struct DDPShocks
	# dimension of aSim is (dimShocks, nPeriods, nFirms)
	# dimension of aInit is (dimShocks, 1, nFirms)
	aInit ::Array{Float64,3}
	aSim ::Array{Float64,3}
end

function drawshocks(p::DDP; nPeriods::Int64 = p.params.nPeriods,
		nFirms::Int64 = p.params.nFirms)

	# dimShocks = size(p.mShocks,1)
	# aInit = randn((dimShocks, 1, nFirms))
	# aSim = randn((dimShocks, nPeriods, nFirms))

	dimShocks = length(p.shockdist)

	aInit = zeros((dimShocks, 1, nFirms))
	for f = 1:nFirms
		aInit[:,1,f] .= rand(p.shockdist)
	end

	aSim = zeros((dimShocks, nPeriods, nFirms))
	for t=1:nPeriods, f=1:nFirms
		aSim[:,t,f] .= rand(p.shockdist)
	end

	return DDPShocks(aInit, aSim)
end

# function drawshocks(p::NewIdeas)
# 	@unpack nFirms, nPeriods, pr_i = p.params
# 	dimShocks = size(p.mShocks,1)
#
# 	aInit_z = randn((1, 1, nFirms))
# 	aInit_i = rand(Float64, (1, 1, nFirms))
# 	aInit = [aInit_z; aInit_i]
#
# 	aSim_z = randn((1, nPeriods, nFirms))
# 	aSim_i = rand(Float64, (1, nPeriods, nFirms))
# 	aSim = [aSim_z; aSim_i]
#
# 	return DDPShocks(aInit, aSim)
# end
