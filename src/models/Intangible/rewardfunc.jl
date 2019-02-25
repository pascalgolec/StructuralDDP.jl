function grossprofits(p::Intangible, vStateVars::Vector{Float64})
	@unpack α, γ = p.params
	(K, N, z) = vStateVars
	return (K^α * N^(1-α))^γ * exp(z)
end

function adjustcosts(p::Intangible,
                     Kprime::Float64,
                     K::Float64,
					 Nprime::Float64,
					 N::Float64)
	@unpack δ_K, δ_N, a_K, a_N = p.params

    # specifiaction follows peters and taylor 2017
    capx = Kprime - (1-δ_K)*K
    intx = Nprime - (1-δ_N)*N

    Ktot = K + N

    adj_K = a_K/2*(capx/Ktot)^2 * Ktot
    adj_N = a_N/2*(intx/Ktot)^2 * Ktot

    return adj_K + adj_N
end

function oibdp(p::Intangible, vStateVars::Vector{Float64}, choice::Vector{Float64})
    @unpack δ_N = p.params
    (K, N, z) = vStateVars
    (Kprime, Nprime) = choice

    intx = Nprime - (1-δ_N)*N
    return grossprofits(p, vStateVars) - adjustcosts(p, Kprime, K, Nprime, N) - intx
end

ebit(p::Intangible, vStateVars::Vector{Float64}, choice::Vector{Float64}) =
    oibdp(p, vStateVars, choice) - p.params.δ_K*vStateVars[1]

function cashflow(p::Intangible, vStateVars::Vector{Float64}, choice::Vector{Float64})
    @unpack δ_K = p.params
    K = vStateVars[1]
    Kprime = choice[1]

    capx = Kprime - (1-δ_K)*K
    return netincome(p, vStateVars, choice) + δ_K * K - capx
end

# for rewardmat prebuild_partial
function rewardfunc(p::Intangible, Output::Float64, K::Float64, Kprime::Float64,
                                                    N::Float64, Nprime::Float64)
    @unpack δ_K, δ_N, τ = p.params

    # specification follows peters and taylor 2017
    capx = Kprime - (1-δ_K)*K
    intx = Nprime - (1-δ_N)*N

    adj = adjustcosts(p, Kprime, K, Nprime, N)

    return (Output - adj - intx)*(1-τ) - capx + δ_K*τ*K
end
