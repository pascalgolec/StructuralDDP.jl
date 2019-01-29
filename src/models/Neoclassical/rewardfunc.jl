####################################
# need this for state-action formulation
####################################
grossprofits(p::NeoClassicalSimple, vStateVars::Vector{Float64}) =
    vStateVars[1]^p.params.α * exp(vStateVars[2])

function adjustcosts(p::NeoClassicalSimple, Kprime::Float64, K::Float64)
    @unpack δ, F, γ = p.params
    capx = Kprime - (1-δ)*K
    action = (Kprime != K)
    return F*K*action + γ/2*(capx/K- δ)^2 * K
end

function oibdp(p::NeoClassicalSimple, vStateVars::Vector{Float64}, Kprime::Float64)
    return grossprofits(p, vStateVars) - adjustcosts(p, Kprime, vStateVars[1])
end

function ebit(p::NeoClassicalSimple, vStateVars::Vector{Float64}, choice)
    return oibdp(p, vStateVars, choice) - p.params.δ*vStateVars[1] # depreciation tax shield
end

netincome(p::DDM, vStateVars::Vector{Float64}, choice) =
                                            ebit(p, vStateVars, choice)*(1-p.params.τ)

function cashflow(p::NeoClassicalSimple, vStateVars::Vector{Float64}, Kprime::Float64)
    @unpack δ, κ = p.params
    K = vStateVars[1]
    capx = Kprime - (1-δ)*K
    return netincome(p, vStateVars, Kprime) + δ*K - (1-κ*(capx<0))*capx
end

dividend(p::DDM, vStateVars::Vector{Float64}, choice) =
    cashflow(p, vStateVars, choice)

rewardfunc(p::DDM, vStateVars::Vector{Float64}, choice) =
    dividend(p, vStateVars, choice)


####################################
# for convexity VFI: also has output and K as input
####################################
function rewardfunc(p::NeoClassicalSimple, Output::Float64, K::Float64, Kprime::Float64)
    @unpack δ, τ, F, γ, κ  = p.params
    capx = Kprime - (1-δ)K
    out = (1-τ)*(Output - F*K - γ/2*(capx/K- δ)^2 * K) - (1-κ*(capx<0))*capx + τ * δ * K
    # dont have condition on F in here!!!!
end

function rewardfuncinaction(p::NeoClassicalSimple, Output::Float64, K::Float64)
    @unpack δ, τ = p.params
    capx = δ*K
    return (1-τ)*Output - capx + τ * δ * K
end
