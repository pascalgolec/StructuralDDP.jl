
####################################################
#### Quantecon state-action representation #########
####################################################

# this function sets up the state-action pairs and feeds them into the QuantEcon solver
function _solve(p::DDM, method::Type{SA},  mTransition, mReward, disp::Bool)

    # ONLY PROGRAMMED FOR SINGLECHOICEVAR SO FAR!

    # state action pair representation
    nStateVars = length(p.tStateVectors)
    tNodes = length.(p.tStateVectors)

    nChoiceVars = length(p.tChoiceVectors)
    nChoices = prod(length.(p.tChoiceVectors))
    # nChoices = length(p.tChoiceVectors[1])

    # tAll = (tNodes..., nChoices)
    # mSA = gridmake(p.tStateVectors..., p.vChoiceVector)

    mSA_coord = gridmake(1:prod(tNodes), 1:nChoices)
    mS_coord = mSA_coord[:,1:1]
    mA_coord = mSA_coord[:,2:2]

    # rewardmatrix gives choices x states, want states x choices
    R = mReward' # is choices x states
    # R_sa = addDim(R[:]) # can convert into state-action repr
    R_sa = R[:]

    # get transition matrix, which elements are allowed?
    Q = mTransition

    # is state-action representation, since the choices don't need to be equal
    # to the states, it is possible that some choices are not admissible
    # this will be the case where probabilities of the transition matrix are not
    # between 0 and 1
    # mark all rows where probabilities are within 0 and 1
    admissible = all((Q.>=0.) .& (Q.<=1.), dims=2)[:]
    # would be more effcient to do this on g rather than Q, huge matrix

    # delete state-action pairs that are not admissible
    mS_coord = mS_coord[admissible]
    mA_coord = mA_coord[admissible]
    R_sa = R_sa[admissible]
    # Q = Q[find(admissible),:]
    Q = Q[(LinearIndices(admissible))[findall(admissible)],:]


    # R_sa = reward
    # Q = transition probability array
    # mS_coord: Action Indices
    # mA_coord: Action Index Pointers
    ddp = DiscreteDP(R_sa, Q, p.params.β, mS_coord, mA_coord)

    results = solve(ddp, VFI)

    if disp
        println(string("number of iterations = ", results.num_iter))
    end

    # results.sigma is a vector of indices of optimal choices
    # mPolFun = p.tChoiceVectors[1][results.sigma]
    vChoiceVector = gridmake(p.tChoiceVectors...)
    mPolFun = vChoiceVector[results.sigma, :]

    # meshPolFun = reshape(mPolFun, tuple(tNodes...))
    meshValFun = reshape(results.v, tuple(tNodes...))

    # return createsolution(p, meshValFun, (meshPolFun,))
    meshPolFun = Vector{Array{Float64}}()
    for k = 1 : nChoiceVars
        push!(meshPolFun, reshape(mPolFun[:,k], tuple(tNodes...)))
    end

    return meshValFun, Tuple(meshPolFun)
end

# function solve(p::TwoChoiceVar, method::Type{SA},  mTransition, mReward, disp::Bool)
#
#     # state action pair representation
#     nStateVars = length(p.tStateVectors)
#     tNodes = length.(p.tStateVectors)
#
#     nChoiceVars = length(p.tChoiceVectors)
#     nChoices = prod(length.(p.tChoiceVectors))
#
#     # tAll = (tNodes..., nChoices)
#     # mSA = gridmake(p.tStateVectors..., p.vChoiceVector)
#
#     mSA_coord = gridmake(1:prod(tNodes), 1:nChoices)
#     mS_coord = mSA_coord[:,1:1]
#     mA_coord = mSA_coord[:,2:2]
#
#     # reward usually gives choices x states, want states x choices
#     R = mReward' # is choices x states
#     # R_sa = addDim(R[:]) # can convert into state-action repr
#     R_sa = R[:]
#
#
#     # get transition matrix, which elements are allowed?
#     Q = transitionmatrix(p, method)
#
#     # mark all rows where probabilities are not within 0 and 1
#     admissible = all((Q.>=0.) .& (Q.<=1.), dims=2)[:]
#     # would be more effcient to do this on g than Q, huge matrix
#
#     # delete s-a pairs that are not admissible
#     mS_coord = mS_coord[admissible]
#     mA_coord = mA_coord[admissible]
#     R_sa = R_sa[admissible]
#     Q = Q[(LinearIndices(admissible))[findall(admissible)], :]
#
#
#     # R_sa = reward
#     # Q = transition probability array
#     # mS_coord: Action Indices
#     # mA_coord: Action Index Pointers
#     ddp = DiscreteDP(R_sa, Q, p.params.β, mS_coord, mA_coord)
#
#     results = solve(ddp, VFI)
#
#     if disp
#         println(string("number of iterations = ", results.num_iter))
#     end
#
#     # size of results.sigma is (nK x nz) x 1
#     # values are the indeces of choices
#     # @show size(results.sigma)
#     # @show results.sigma[1:10]
#
#     vChoiceVector = gridmake(p.tChoiceVectors...)
#     mPolFun = vChoiceVector[results.sigma, :]
#
#     meshPolFun1 = reshape(mPolFun[:,1], tuple(tNodes...))
#     meshPolFun2 = reshape(mPolFun[:,2], tuple(tNodes...))
#     meshValFun = reshape(results.v, tuple(tNodes...))
#
#     return createsolution(p, meshValFun, (meshPolFun1, meshPolFun2))
#
# end
