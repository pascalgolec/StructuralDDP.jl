# construct matrices of indirect utility and initial choice depending on exogenous variables
# so in the neoclassical model, only have a grid with K0 as a function of z0


#######################
# For state-action representation
#######################
# to do:
# - only works if choice vector is equal to state vector
# - only works if first state is choice state?
function initialendogstatevars(p::SingleChoiceVar, meshValFun)

    tNodes = length.(p.tStateVectors)
    exogtNodes = tNodes[.!p.bEndogStateVars]

    mV0 = zeros(exogtNodes)
    mChoice0 = zeros(exogtNodes)

    # vChoice = p.vChoiceVector
    vInitialState = p.tStateVectors[1]

    for iexog in CartesianIndices(exogtNodes)

        interimvalue = -Inf
        interimChoice = 0

        for iChoice = 1:length(vInitialState)

                # fundsneeded =  (1+p.C0)*vK[iK] - vLev[iLev]*vK[iK]
                # fincf = dividends(-fundsneeded, p) # fincf is negative
                # value = meshValFun[iChoice,iexog.I...] - (1+p.C0)* vChoice[iChoice]

                value = initializationproblem(p, meshValFun[iChoice,iexog.I...], vInitialState[iChoice])

                if value > interimvalue
                    interimvalue = value
                    interimChoice = vInitialState[iChoice]
                end

        end

        mV0[iexog.I...] = interimvalue
        mChoice0[iexog.I...] = interimChoice

    end

    return mV0, mChoice0

end

# function initialendogstatevars(p::ManyChoiceVars, meshValFun)
#
#     tNodes = length.(p.tStateVectors)
#     exogtNodes = tNodes[.!p.bEndogStateVars]
#
#     mV0 = zeros(exogtNodes)
#     mChoiceOne0 = zeros(exogtNodes)
#     mChoiceTwo0 = zeros(exogtNodes)
#
#     vChoiceOne = p.tChoiceVectors[1]
#     vChoiceTwo = p.tChoiceVectors[2]
#
#     for iexog in CartesianIndices(exogtNodes)
#
#         interimvalue = -Inf
#         interimChoiceOne = 0
#         interimChoiceTwo = 0
#
#         for iChoiceOne = 1:length(vChoiceOne), iChoiceTwo = 1:length(vChoiceTwo)
#
#                 # fundsneeded =  (1+p.C0)*vK[iK] - vLev[iLev]*vK[iK]
#                 # fincf = dividends(-fundsneeded, p) # fincf is negative
#                 # value = meshValFun[iChoice,iexog.I...] - (1+p.C0)* vChoice[iChoice]
#
#                 value = initializationproblem(p, meshValFun[iChoiceOne, iChoiceTwo, iexog.I...],
#                             vChoiceOne[iChoiceOne], vChoiceTwo[iChoiceTwo])
#
#                 if value > interimvalue
#                     interimvalue = value
#                     interimChoiceOne = vChoiceOne[iChoiceOne]
#                     interimChoiceTwo = vChoiceTwo[iChoiceTwo]
#                 end
#
#         end
#
#         mV0[iexog.I...] = interimvalue
#         mChoiceOne0[iexog.I...] = interimChoiceOne
#         mChoiceTwo0[iexog.I...] = interimChoiceTwo
#
#     end
#
#     return mV0, mChoiceOne0, mChoiceTwo0
#
# end
