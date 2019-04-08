

function getquadrature(d::ContinuousUnivariateDistribution, numquadnodes::Vector{Int})
	vShocks, vWeights = qnwdist(d, numquadnodes[1])
    mShocks = Array(vShocks')
    return mShocks, vWeights
end
function getquadrature(d::Normal, numquadnodes::Vector{Int})
    vShocks, vWeights = qnwnorm(numquadnodes[1], mean(d), std(d)^2)
    mShocks = Array(vShocks')
    return mShocks, vWeights
end
function getquadrature(d::Uniform, numquadnodes::Vector{Int})
    vShocks, vWeights = qnwunif(numquadnodes[1], params(d)...)
    mShocks = Array(vShocks')
    return mShocks, vWeights
end
function getquadrature(d::LogNormal, numquadnodes::Vector{Int})
    vShocks, vWeights = qnwlogn(numquadnodes[1], meanlogx(d), varlogx(d))
    mShocks = Array(vShocks')
    return mShocks, vWeights
end

getquadrature(d::MvNormal, numquadnodes::Vector{Int}) =
	mShocks, vWeights = qnwnorm(numquadnodes, Vector(d.μ), Matrix(d.Σ))

# split into intdims
transitionmatrix(p::DDM;
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ))) =
    _transitionmatrix(p, p.intdim, numquadnodes)


"""Generate a nExogStates x nExogStates transition matrix."""
_transitionmatrix(p::DDM, method::Type{Separable_ExogStates}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
	getnonchoicevars(p.tStateVectors, p.tChoiceVectors),
	getnonchoicevars(p.tStateVectors, p.tChoiceVectors),
    p.shockdist, numquadnodes)

"""Generate a nStates x nExogStates transition matrix."""
_transitionmatrix(p::DDM, method::Type{Separable_States}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
    p.tStateVectors,
	getnonchoicevars(p.tStateVectors, p.tChoiceVectors),
   	p.shockdist, numquadnodes)

"""Generate a (nChoices*nStates) x nExogStates transition matrix."""
_transitionmatrix(p::DDM, method::Type{Separable}, numquadnodes::Vector{Int}) =
    _transitionmatrix(p.transfunc,
	    p.tStateVectors, getchoicevars(p.tStateVectors, p.tChoiceVectors),
		getnonchoicevars(p.tStateVectors, p.tChoiceVectors),
	    p.shockdist, numquadnodes)

"""Generate a (nChoices*nStates) x nStates transition matrix."""
_transitionmatrix(p::DDM, method::Type{All}, numquadnodes::Vector{Int}) =
	_transitionmatrix(p.transfunc,
		p.tStateVectors, p.tChoiceVectors,
		p.tStateVectors,
		p.shockdist, numquadnodes)


# for choices and states as input
function _transitionmatrix(transfunc::Function,
    tInputVectorsStates::NTuple{N1, Vector{T}},
	tInputVectorsChoices::NTuple{N2, Vector{T}},
	tOutputVectors,
    shockdist, numquadnodes) where {N1, N2, T}

    mShocks, vWeights = getquadrature(shockdist, numquadnodes)

    nOutputStates = prod(length.(tOutputVectors))
    dimOutputStates = length(tOutputVectors)

    # nChoices = prod(length.(tChoiceVectors))
    # nChoices = prod(length.(tChoiceVectors))
    nInputStates = prod(length.(tInputVectorsStates)) * prod(length.(tInputVectorsChoices))

    g = zeros(dimOutputStates, nInputStates)
    mG = spzeros(nInputStates, nOutputStates)
    basisOutputStates = Basis(SplineParams.(tOutputVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

	# choices change faster than states
	iterator = Iterators.product(Iterators.product(tInputVectorsChoices...), Iterators.product(tInputVectorsStates...))
	# iterator = Iterators.product(Iterators.product(tInputVectorsStates...), Iterators.product(tInputVectorsChoices...))

	# loop over all shock combinations
    for i = 1 : length(vWeights)
        # loop over all inputstates
        for (j, (choices, states)) in enumerate(iterator)
        # for (j, (states, choices)) in enumerate(iterator)
			# @show j
			# @show choices
			# @show states
            g[:,j] .= transfunc(states, choices,  mShocks[:,i])
			# @show g[:,j]
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

	mShocks, vWeights = getquadrature(shockdist, numquadnodes)

	nInputStates = prod(length.(tInputVectors))
    dimOutputStates = length(tOutputVectors)
    nOutputStates = prod(length.(tOutputVectors))

    g = zeros(dimOutputStates, nInputStates)
    mG = spzeros(nInputStates, nOutputStates)
    basisOutputStates = Basis(SplineParams.(tOutputVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

	iterator = Iterators.product(tInputVectors...)

    # loop over all shock combinations
    for i = 1 : length(vWeights)

        # loop over all inputstates
        for (j, states) in enumerate(iterator)
            g[:,j] .= transfunc(states, mShocks[:,i])
        end

        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + vWeights[i] * PhiTemp
    end
    return mG::SparseMatrixCSC{Float64,Int64}

end
