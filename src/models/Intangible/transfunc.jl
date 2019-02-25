function transfunc(p::Intangible, method::Type{separable}, vExogState, vShocks::Vector{Float64})
    @unpack ρ, σ  = p.params

    z = vExogState[1]
    zprime  = ρ*z + σ * vShocks[1];
    return  inbounds(zprime, p.tStateVectors[3][1], p.tStateVectors[3][end])
end

transfunc(p::Intangible, method::Type{intermediate}, vState, vShocks) =
    transfunc(p, separable, vState[3], vShocks)

transfunc(p::Intangible, method::Type{SA}, vState, vChoice, vShock) =
    [vChoice[1], vChoice[2], transfunc(p, intermediate, vState, vShock)]
