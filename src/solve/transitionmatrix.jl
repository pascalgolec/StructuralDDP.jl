

getquadrature(d::ContinuousUnivariateDistribution, numquadnodes::Vector{Int}) =
	qnwdist(d, numquadnodes[1])

getquadrature(d::Normal, numquadnodes::Vector{Int}) =
    qnwnorm(numquadnodes[1], mean(d), std(d)^2)

getquadrature(d::Uniform, numquadnodes::Vector{Int}) =
    qnwunif(numquadnodes[1], params(d)...)

getquadrature(d::LogNormal, numquadnodes::Vector{Int}) =
    qnwlogn(numquadnodes[1], meanlogx(d), varlogx(d))

getquadrature(d::MvNormal, numquadnodes::Vector{Int}) =
	qnwnorm(numquadnodes, Vector(d.μ), Matrix(d.Σ))


transitionmatrix(p::DDP;
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ))) =
    _transitionmatrix(p, p.intdim, numquadnodes)


"""Generate a nExogStates x nExogStates transition matrix."""
_transitionmatrix(p::DDP, method::Type{Separable_ExogStates}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
	getnonchoicevars(p),
	getnonchoicevars(p),
    p.shockdist, numquadnodes)

"""Generate a nStates x nExogStates transition matrix."""
_transitionmatrix(p::DDP, method::Type{Separable_States}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
    p.tStateVectors,
	getnonchoicevars(p),
   	p.shockdist, numquadnodes)

"""Generate a (nChoices*nStates) x nExogStates transition matrix."""
_transitionmatrix(p::DDP, method::Type{Separable}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
	    p.tStateVectors, getchoicevars(p),
		getnonchoicevars(p),
	    p.shockdist, numquadnodes)

"""Generate a (nChoices*nStates) x nStates transition matrix."""
_transitionmatrix(p::DDP, method::Type{All}, numquadnodes::Vector{Int}) =
	_transitionmatrix(p.transfunc,
		p.tStateVectors, p.tChoiceVectors,
		p.tStateVectors,
		p.shockdist, numquadnodes)

# for choices and states as input
function _transitionmatrix(transfunc::Function,
    tInputVectorsStates::NTuple{dimStates, Vector{T1}},
	tInputVectorsChoices::NTuple{dimChoices, Vector{T2}},
	tOutputVectors,
    shockdist, numquadnodes) where {dimStates, dimChoices, T1, T2}

    Shocks, vWeights = getquadrature(shockdist, numquadnodes)
	dimshocks = length(shockdist)

    nOutputStates = prod(length.(tOutputVectors))
    dimOutputStates = length(tOutputVectors)

    nInputStates = prod(length.(tInputVectorsStates)) * prod(length.(tInputVectorsChoices))

    g = zeros(dimOutputStates, nInputStates)
    mG = spzeros(nInputStates, nOutputStates)
    basisOutputStates = Basis(SplineParams.(tOutputVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

	iterator_choices = getiterator(tInputVectorsChoices)

	# choices change faster than states
	iterator = Iterators.product(iterator_choices, Iterators.product(tInputVectorsStates...))

	lowerbounds = minimum.(tOutputVectors)
	upperbounds = maximum.(tOutputVectors)

	# loop over all shock combinations
    for i = 1 : length(vWeights)
        # loop over all inputstates
        for (j, (choices, states)) in enumerate(iterator)
			if dimshocks > 1
            	# g[:,j] .= inbounds.(transfunc(states, choices,  Shocks[:,i]),
				# 	lowerbounds, upperbounds)
            	g[:,j] .= transfunc(states, choices,  Shocks[:,i])
			else
				# g[:,j] .= inbounds.(transfunc(states, choices,  Shocks[i]),
				# 	lowerbounds, upperbounds)
				g[:,j] .= transfunc(states, choices,  Shocks[i])
			end
        end
        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + vWeights[i] * PhiTemp
    end
    # is G fully sparse?
    # G[G < 1E-10] = 0
    # dropzeros!(G)

    # VERY IMPORTANT TO SPECIFY TYPE HERE, BECAUSE BASISMATRICES DOES NOT GIVE FIXED OUTPUT
    #  if gettransitionmatrix turns out to be a bottleneck, then I need to dig deeper into this
    return mG::SparseMatrixCSC{Float64,Int64}
end #getG

# for states as input
function _transitionmatrix(transfunc::Function,
    tInputVectors::NTuple{N, Vector{T}}, tOutputVectors,
    shockdist, numquadnodes) where {N,T}

	Shocks, vWeights = getquadrature(shockdist, numquadnodes)
	dimshocks = length(shockdist)

	nInputStates = prod(length.(tInputVectors))
    dimOutputStates = length(tOutputVectors)
    nOutputStates = prod(length.(tOutputVectors))

    g = zeros(dimOutputStates, nInputStates)
    mG = spzeros(nInputStates, nOutputStates)
    basisOutputStates = Basis(SplineParams.(tOutputVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

	iterator = getiterator(tInputVectors)

	lowerbounds = minimum.(tOutputVectors)
	upperbounds = maximum.(tOutputVectors)

    # loop over all shock combinations
    for i = 1 : length(vWeights)
        # loop over all inputstates
        for (j, states) in enumerate(iterator)
			if dimshocks > 1
            	# g[:,j] .= inbounds.(transfunc(states, Shocks[:,i]),
				# 	lowerbounds, upperbounds)
            	g[:,j] .= transfunc(states, Shocks[:,i])
        	else
				# g[:,j] .= inbounds.(transfunc(states, Shocks[i]),
				# 	lowerbounds, upperbounds)
				g[:,j] .= transfunc(states, Shocks[i])
			end
        end
        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + vWeights[i] * PhiTemp
    end
    return mG::SparseMatrixCSC{Float64,Int64}

end
