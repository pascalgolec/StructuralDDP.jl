
function initializationproblem(p::NeoClassicalSimple, value::Float64, K::Float64)
    @unpack β, C0 = p.params
    r = (1-β)/β
    return value - (1 + r + C0) * K
end

function initialize(p::NeoClassicalSimple, dShock::Array{Float64, 1}, itp_K0)
	@unpack σ, ρ = p.params
    dShock = dShock[1]
    z0 = dShock * sqrt(σ^2 / (1-ρ^2))
	# z0 = inbounds(z0, p.vMin[2], p.vMax[2])
	z0 = inbounds(z0, p.tStateVectors[2][1], p.tStateVectors[2][end])

    K0 = itp_K0(z0)
	K0 = inbounds(K0, p.tStateVectors[1][1], p.tStateVectors[1][end])
	return [K0, z0]
end
