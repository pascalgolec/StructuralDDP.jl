
transitionmatrix(p::DDM; intdim::Symbol=:SA,
    numquadnodes::Vector{Int} = 5*ones(Int64, length(p.shockdist.μ))) =
    _transitionmatrix(p, eval(intdim), numquadnodes)

# _transitionmatrix(p::DDM, method::Type{T}) where T<:DDMIntDim =
#     _transitionmatrix(method, p.transfunc, p.tStateVectors, p.tChoiceVectors,
#     p.bEndogStateVars,
#     p.vWeights, p.mShocks)

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



# expand structure, only use necessary inputs
_transitionmatrix(p::DDM, method::Type{separable}, numquadnodes::Vector{Int}) =
    _transitionmatrix(method, p.transfunc, p.tStateVectors[.!p.bEndogStateVars],
    p.shockdist, numquadnodes)
function _transitionmatrix(method::Type{separable}, transfunc::Function,
    tStochStateVectors::NTuple{N,Vector{Float64}},
    shockdist, numquadnodes) where N
    # Calculates probability transition matrix, depending on integration dimension
    # - uses gaussian quadrature
    # - uses linear basis functions
    # - size of matrix is nInputStates x nOutputStates.
    #       - InputStates together with shocks affect the transition of OutputStates
    #       - so it is States_today x States_tomorrow, same dimension as QuantEcon.tauchen()

    mShocks, vWeights = getquadrature(shockdist, numquadnodes)

    dimStochStates = length(tStochStateVectors)
    nStochStates = prod(length.(tStochStateVectors))

    # all combinations of exogenous variables
    mStochStates = gridmake(tStochStateVectors...)
    if length(size(mStochStates)) == 1 # need to add dimension
        mStochStates = addDim(mStochStates)
    end
    mStochStates = mStochStates'

    g = zeros(dimStochStates, nStochStates)
    mG = spzeros(nStochStates, nStochStates)
    basisOutputStates = Basis(SplineParams.(tStochStateVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

    # loop over all shock combinations
    for i = 1 : length(vWeights)

        # loop over all inputstates
        for j = 1 : nStochStates
            g[:,j] .= transfunc(method, mStochStates[:, j], mShocks[:,i])
            # can speed up slightly if only have one shock
            # g[:,j] .= transfunc(method, mStochStates[:, j], mShocks[1,i])
        end

        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + vWeights[i] * PhiTemp
    end
    return mG::SparseMatrixCSC{Float64,Int64}
end

_transitionmatrix(p::DDM, method::Type{intermediate}, numquadnodes::Vector{Int}) =
    _transitionmatrix(method, p.transfunc,
    p.tStateVectors, p.tStateVectors[.!p.bEndogStateVars],
   p.shockdist, numquadnodes)
function _transitionmatrix(method::Type{intermediate}, transfunc::Function,
    tStateVectors, tStochStateVectors,
    shockdist, numquadnodes)
    # mG is nStates x nStochStates

    mShocks, vWeights = getquadrature(shockdist, numquadnodes)

    nStates = prod(length.(tStateVectors))
    mStates = gridmake(tStateVectors...)
    if length(size(mStates)) == 1 # need to add dimension
        mStates = addDim(mStates)
    end
    mStates = mStates'

    dimStochStates = length(tStochStateVectors)
    nStochStates = prod(length.(tStochStateVectors))

    g = zeros(dimStochStates, nStates)
    mG = spzeros(nStates, nStochStates)
    basisOutputStates = Basis(SplineParams.(tStochStateVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

    # loop over all shock combinations
    for i = 1 : length(vWeights)
        # loop over all inputstates
        for j = 1 : nStates
            g[:,j] .= transfunc(method, mStates[:, j], mShocks[:,i])
            # g[:,j] .= transfunc(method, mStates[:, j], mShocks[1,i])
        end

        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + vWeights[i] * PhiTemp
    end
    return mG::SparseMatrixCSC{Float64,Int64}
end #getG

_transitionmatrix(p::DDM, method::Type{SA}, numquadnodes::Vector{Int}) =
    _transitionmatrix(method, p.transfunc,
    p.tStateVectors, p.tChoiceVectors,
    p.shockdist, numquadnodes)
function _transitionmatrix(method::Type{SA}, transfunc::Function,
    tStateVectors, tChoiceVectors,
    shockdist, numquadnodes)
    # mG is nChoices*nStates x nStates

    mShocks, vWeights = getquadrature(shockdist, numquadnodes)

    nStates = prod(length.(tStateVectors))
    dimStates = length(tStateVectors)
    nChoices = prod(length.(tChoiceVectors))
    nInputStates = nStates * nChoices

    mInputStates = gridmake(tStateVectors..., tChoiceVectors...)
    mInputStates_states = mInputStates[:, 1:dimStates]
    vInputStates_choices = mInputStates[:, dimStates+1 : end]

    g = zeros(dimStates, nInputStates)
    mG = spzeros(nInputStates, nStates)
    basisOutputStates = Basis(SplineParams.(tStateVectors,0,1))
    PhiTemp = spzeros(size(mG)...)

    # loop over all shock combinations
    for i = 1 : length(vWeights)
        # loop over all inputstates
        for j = 1 : nInputStates

            g[:,j] .= transfunc(method, mInputStates_states[j, :],
                vInputStates_choices[j, :],  mShocks[:,i])

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
