# TRANSFUNC
# - separable: exogenous state variables as input
# - intermediate: all state variables as input
# - SA: all state variables and choice variables as input

#############
# NEOCLASSICAL
#############
function transfunc(p::NeoClassicalSimple, method::Type{separable}, vExogState, vShocks)
    @unpack ρ , σ = p.params
    z = vExogState[1]
    zprime  = ρ*z + σ * vShocks[1];
    return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
end
#
# transfunc(p::NeoClassicalSimple, method::Type{intermediate}, vState, vShocks) =
#     transfunc(p, separable, vState[2], vShocks)
#
# transfunc(p::NeoClassicalSimple, method::Type{SA}, vState, vChoice, vShock) =
#     [vChoice..., transfunc(p, intermediate, vState, vShock)...]

# function transfunc(p::NeoClassicalSimple,
#     vExogState::Vector{Float64}, vChoice::Nothing, vShock::Vector{Float64})
#     @unpack ρ , σ = p.params
#     z = vState[2]
#     zprime  = ρ*z + σ * vShock[1];
#     return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
# end
