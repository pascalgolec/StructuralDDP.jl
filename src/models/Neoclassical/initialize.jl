

function initialize(p::NeoClassicalSimple, dShock::Array{Float64, 1}, itp_K0)
	@unpack σ, ρ = p.params
    dShock = dShock[1]
    z0 = dShock* sqrt(σ^2 / (1-ρ^2))
	z0 = inbounds(z0, p.vMin[2], p.vMax[2])

    K0 = itp_K0(z0)
	K0 = inbounds(K0, p.vMin[1], p.vMax[1])
	return [K0, z0]
end
