
transitionmatrix(p::DDM; intdim::Symbol=:SA) = transitionmatrix(p, eval(intdim))

# if user only supplies separable, can still do the others with the help of these
transfunc(p::DDM, method::Type{intermediate}, vState, vShocks) =
    transfunc(p, separable, vState[.!p.bEndogStateVars], vShocks)
transfunc(p::DDM, method::Type{SA}, vState, vChoice, vShock) =
    [vChoice..., transfunc(p, intermediate, vState, vShock)...]

# transfunc(p::DDM, method::Type{separable}, vState, vShock) =
#     transfunc(p::DDM, vState, nothing, vShock)
# transfunc(p::DDM, method::Type{intermediate}, vState, vShock) =
#     transfunc(p::DDM, vState, nothing, vShock)

function transitionmatrix(p::DDM, method::Type{T}) where T<:DDMIntDim
    # Calculates probability transition matrix, depending on integration dimension
    # - uses gaussian quadrature
    # - uses linear basis functions
    # - size of matrix is nInputStates x nOutputStates.
    #       - InputStates together with shocks affect the transition of OutputStates
    #       - so it is States_today x States_tomorrow, same dimension as QuantEcon.tauchen()

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

            # if method == SA
            #     g[:,j] .= transfunc(p, mInputStates_states[j, :],
            #         vInputStates_choices[j, :],  p.mShocks[:,i])
            # else
            #     g[:,j] .= transfunc(p, mInputStates[:, j], p.mShocks[:,i])
            # end
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
