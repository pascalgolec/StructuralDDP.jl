# Note: shocks should be independent of parameters, because in SMM will loop over
# parameters, keeping the shocks fixed

struct DDPShocks
	# dimension of aSim is (dimShocks, nPeriods, nFirms)
	# dimension of aInit is (dimShocks, 1, nFirms)
	aInit ::Array{Float64,3}
	aSim ::Array{Float64,3}
end

# drawshocks(p::DDM) = drawshocks(p, p.params.nPeriods, p.params.nFirms)

function drawshocks(p::DDM; nPeriods::Int64 = p.params.nPeriods,
		nFirms::Int64 = p.params.nFirms)
	# @unpack nFirms, nPeriods = p.params
	dimShocks = size(p.mShocks,1)
	aInit = randn((dimShocks, 1, nFirms))
	aSim = randn((dimShocks, nPeriods, nFirms))

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
