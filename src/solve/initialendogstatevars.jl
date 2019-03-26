# construct matrices of indirect utility and initial choice depending on exogenous variables
# so in the neoclassical model, only have a grid with K0 as a function of z0


#######################
# For state-action representation
#######################
# to do:
# - only works if choice vector is equal to state vector
# - only works if first state is choice state?
# note: could be more elegant to have a SIngleChoiceVar and a TwoChoiceVar subtype of DDM
# to allocate to method.. will depend on how user inputs stuff
function initialendogstatevars(p::DDM, meshValFun::Array{Float64})
    if typeof(p.tChoiceVectors) == Tuple{Vector{Float64}}
        return meshValFunZero, tmeshPolFunZero = initialendogstatevars1(
            p.initializationproblem, meshValFun, p.tStateVectors, p.bEndogStateVars)
    elseif typeof(p.tChoiceVectors) == Tuple{Vector{Float64},Vector{Float64}}
        return meshValFunZero, tmeshPolFunZero = initialendogstatevars2(
            p.initializationproblem, meshValFun, p.tStateVectors, p.bEndogStateVars)
    end
end

function initialendogstatevars1(initializationproblem::Function, meshValFun::Array{Float64},
    tStateVectors, bEndogStateVars)

    tNodes = length.(tStateVectors)
    exogtNodes = tNodes[.!bEndogStateVars]

    mV0 = zeros(exogtNodes)
    mChoice0 = zeros(exogtNodes)

    # vActions = tStateVectors[1]
    vActions::Vector{Real} = tStateVectors[bEndogStateVars][1]

    for iexog in CartesianIndices(exogtNodes)

        interimvalue = -Inf
        interimChoice = 0

        for iChoice = 1:length(vActions)

                value = initializationproblem(meshValFun[iChoice,iexog.I...], vActions[iChoice])

                if value > interimvalue
                    interimvalue = value
                    interimChoice = vActions[iChoice]
                end

        end

        mV0[iexog.I...] = interimvalue
        mChoice0[iexog.I...] = interimChoice

    end

    return mV0, (mChoice0,)

end

function initialendogstatevars2(initializationproblem::Function, meshValFun::Array{Float64},
    tStateVectors, bEndogStateVars)

    tNodes = length.(tStateVectors)
    exogtNodes = tNodes[.!bEndogStateVars]

    mV0 = zeros(exogtNodes)
    mChoiceOne0 = zeros(exogtNodes)
    mChoiceTwo0 = zeros(exogtNodes)

    # vActionsOne = tChoiceVectors[1]
    # vActionsTwo = tChoiceVectors[2]
    vActionsOne, vActionsTwo = tStateVectors[bEndogStateVars]


    for iexog in CartesianIndices(exogtNodes)

        interimvalue = -Inf
        interimChoiceOne = 0
        interimChoiceTwo = 0

        for iChoiceOne = 1:length(vActionsOne), iChoiceTwo = 1:length(vActionsTwo)

                value = initializationproblem(meshValFun[iChoiceOne, iChoiceTwo, iexog.I...],
                            vActionsOne[iChoiceOne], vActionsTwo[iChoiceTwo])

                if value > interimvalue
                    interimvalue = value
                    interimChoiceOne = vActionsOne[iChoiceOne]
                    interimChoiceTwo = vActionsTwo[iChoiceTwo]
                end

        end

        mV0[iexog.I...] = interimvalue
        mChoiceOne0[iexog.I...] = interimChoiceOne
        mChoiceTwo0[iexog.I...] = interimChoiceTwo

    end

    return mV0, (mChoiceOne0, mChoiceTwo0)

end
