
transitionmatrix(p::DDM) = transitionmatrix(p, eval(p.params.intdim))
# transitionmatrix(p::DDM, method::Symbol) = transitionmatrix(p, eval(method))

function transitionmatrix(p::DDM, method::Type{T}) where T<:DDMIntDim
    # Calculates probability transition matrix G of stoch. vars
    # - use gaussian quadrature to get G, not Tauchen method
    # - we use linear basis functions to get the transition matrix

    # size of mTransition will be (nInputStates x nOutputStates)
    # InputStates together with shocks affect the transition of OutputStates

    nStates = prod(length.(p.tStateVectors))
    dimStates = length(p.tStateVectors)

    # if typeof(p) <: SingleChoiceVar
    #     nChoices = length(p.vChoiceVector)
    # else
    #     nChoices = prod(length.(p.tChoiceVectors))
    # end

    nChoices = prod(length.(p.tChoiceVectors))

    bStochStateVars = .!p.bEndogStateVars
    dimStochStates = sum(.!p.bEndogStateVars)
    tStochStateVars = p.tStateVectors[.!p.bEndogStateVars]
    nStochStates = prod(length.(tStochStateVars))

    if method == separable # mT is nStochStates x nStochStates
        nInputStates = nStochStates
        mInputStates = gridmake(tStochStateVars...)
        basisOutputStates = Basis(SplineParams.(tStochStateVars,0,1))
        g = zeros(dimStochStates, nInputStates)
        mG = spzeros(nInputStates, nStochStates)
    elseif method == intermediate # mT is nStates x nStochStates
        nInputStates = nStates
        mInputStates = gridmake(p.tStateVectors...)
        basisOutputStates = Basis(SplineParams.(tStochStateVars,0,1))
        g = zeros(dimStochStates, nInputStates)
        mG = spzeros(nInputStates, nStochStates)
    elseif method == SA # mT is nChoices*nStates x nStates
        nInputStates = nStates * nChoices
        # if typeof(p) <: SingleChoiceVar
        #     mInputStates = gridmake(p.tStateVectors..., p.vChoiceVector)
        # elseif typeof(p) <: ManyChoiceVars
        #     mInputStates = gridmake(p.tStateVectors..., p.tChoiceVectors...)
        # end
        mInputStates = gridmake(p.tStateVectors..., p.tChoiceVectors...)

        mInputStates_states = mInputStates[:, 1:dimStates]
        vInputStates_choices = mInputStates[:, dimStates+1 : end]
        basisOutputStates = Basis(SplineParams.(p.tStateVectors,0,1))
        g = zeros(dimStates, nInputStates)
        mG = spzeros(nInputStates, nStates)
    end

    if length(size(mInputStates)) == 1 # need to add dimension
        mInputStates = addDim(mInputStates)
    end
    mInputStates = mInputStates'

    PhiTemp = spzeros(size(mG)...)

    # loop over all shock combinations
    for i = 1 : length(p.vWeights)

        # loop over all inputstates
        for j = 1 : nInputStates
            if method == SA
                g[:,j] .= transfunc(p, method, mInputStates_states[j, :],
                    vInputStates_choices[j, :],  p.mShocks[:,i])
            else
                g[:,j] .= transfunc(p, method, mInputStates[:, j], p.mShocks[:,i])
            end
        end

        PhiTemp .= BasisMatrix(basisOutputStates, Expanded(), g', 0).vals[1]
        mG .= mG + p.vWeights[i] * PhiTemp

    end

    # is G fully sparse?
    # G[G < 1E-10] = 0
    # dropzeros!(G)

    # VERY IMPORTANT TO SPECIFY TYPE HERE, BECAUSE BASISMATRICES DOES NOT GIVE FIXED OUTPUT
    #  if gettransitionmatrix turns out to be a bottleneck, then I need to dig deeper into this
    return mG::SparseMatrixCSC{Float64,Int64}
end #getG
