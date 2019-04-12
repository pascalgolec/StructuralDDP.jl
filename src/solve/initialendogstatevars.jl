# construct matrices of indirect utility and initial choice depending on exogenous variables
# so in the neoclassical model, only have a grid with K0 as a function of z0


#######################
# For state-action representation
#######################
# to do:
# - only works if choice vector is equal to state vector
# - only works if first state is choice state?
# note: could be more elegant to have a SingleChoiceVar and a TwoChoiceVar subtype of DDM
# to allocate to method.. will depend on how user inputs stuff

initialendogstatevars(p::DDP, meshValFun::Array{Float64}) =
    _initialendogstatevars(
        p.initializationproblem, meshValFun,
        p.tStateVectors, getchoicevarszero(p), getnonchoicevarszero(p))

function _initialendogstatevars(initializationproblem::Function, meshValFun::Array{Float64},
    tStateVectors,
    tChoiceVectors::Tuple{Vector{Float64}}, tExogStateVectors)

    exogtNodes = length.(tExogStateVectors)
    mV0 = zeros(exogtNodes)
    mChoice0 = zeros(exogtNodes)

    vActions::Vector{Real} = tChoiceVectors[1]

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

function _initialendogstatevars(initializationproblem::Function, meshValFun::Array{Float64},
    tStateVectors, tChoiceVectors::NTuple{2,Vector{Float64}}, tExogStateVectors)

    exogtNodes = length.(tExogStateVectors)
    mV0 = zeros(exogtNodes)
    mChoiceOne0 = zeros(exogtNodes)
    mChoiceTwo0 = zeros(exogtNodes)

    vActionsOne, vActionsTwo = tChoiceVectors

    for iexog in CartesianIndices(exogtNodes)

        interimvalue = -Inf
        interimChoiceOne = 0
        interimChoiceTwo = 0

        for iChoiceOne = 1:length(vActionsOne), iChoiceTwo = 1:length(vActionsTwo)

                value = initializationproblem(meshValFun[iChoiceOne, iChoiceTwo, iexog.I...],
                            (vActionsOne[iChoiceOne], vActionsTwo[iChoiceTwo]))

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
