

function initializationproblem(p::Intangible, value::Float64,
                                    choiceone::Float64, choicetwo::Float64)
    @unpack C0 = p.params
    K = choiceone
    N = choicetwo
    # want to choose K to maximize V[K,N,z] - (1+p.C0)*K - (1+p.C0)*N
    return value - (1+C0)*K - (1+C0)*N
end

function initialize(p::Intangible, vShocks::Array{Float64, 1}, itp_K0, itp_N0)
	@unpack σ, ρ = p.params

	dShock = vShocks[1]
    z0 = dShock* sqrt(σ^2 / (1-ρ^2))
	z0 = inbounds(z0, p.tStateVectors[3][1], p.tStateVectors[3][end])

	K0 = itp_K0(z0)
	K0 = inbounds(K0, p.tStateVectors[1][1], p.tStateVectors[1][end])

	N0 = itp_N0(z0)
	N0 = inbounds(N0, p.tStateVectors[2][1], p.tStateVectors[2][end])

	return [K0, N0, z0]
end
