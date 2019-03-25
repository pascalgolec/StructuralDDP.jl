@with_kw struct NeoClassicalSimpleParams{R<:Real, I<:Int64} <: ModelParams
    α::R = 0.67
    β::R = 0.9
    ρ::R = 0.6
    σ::R = 0.3
    δ::R = 0.15
    γ::R = 2.0
    F::R = 0.01
    C0::R = 0.6
    π::R = 0.
    τ::R = 0.35

    κ::R = 0.
    nK::I = 100
    nz::I = 20

    nShocks::I = 3
    nPeriods::I = 120
    nFirms::I = 1000

    # rewardmat::Symbol = :nobuild
    # intdim::Symbol = :separable
    # monotonicity::B = true
    # concavity::B = true
end

function NeoClassicalSimple(; kwargs...)

    # all userprovided parameters go in here
    params = NeoClassicalSimpleParams{Float64, Int64}(; kwargs...)

    # need to unpack to get default values if not provided
    @unpack_NeoClassicalSimpleParams params

    # nNodes = collect((nK, nz))

    # calculate Amin, Amax
    stdz = sqrt(σ^2/(1-ρ^2))
    minz = -3*stdz
    maxz =  3*stdz
    vz = collect(LinRange(minz, maxz, nz))

    function steadystate(z)
        # check https://www3.nd.edu/~esims1/investment_notes.pdf

        function f2!(F,x)
            (V_K, K_log) = x

            F[1] = - 1 - γ * 0 * (1-τ) + V_K
            F[2] = V_K - β*( exp(z)*α*exp(K_log)^(α-1)*(1-τ) + τ*δ + (1-δ)*V_K)
        end

        return nlsolve(f2!, [1.1; 1.])
    end

    K_ss = exp(steadystate(stdz^2/2).zero[2])
    minK_log = steadystate(-2*stdz).zero[2]
    maxK_log = steadystate(3*stdz).zero[2]
    vK   = exp.(collect(LinRange(minK_log, maxK_log, nK)))

    # is incorrect
    # K_ss_analytical(z) = (α * (1-τ) * exp(z)/ ((1-β)/β*(1+γ*δ) + δ)) ^ (1/(1-α))
    # @show log(K_ss_analytical(minz))
    # @show log(K_ss_analytical(maxz))

    # vMin = [exp(minK_log), minz]
    # vMax = [exp(maxK_log), maxz]

    tStateVectors = (vK, vz) # tuple of basis vectors
    # rChoiceStateVars = 1:1

    tChoiceVectors = (vK,) # important to add the comma, otherwise not a tuple

    bEndogStateVars = [true, false]

    # shocks
    vShocks, vWeights = qnwnorm(nShocks,0,1)
    mShocks = vShocks'

    # NeoClassicalSimple(params, tStateVectors, tChoiceVectors,
    #     bEndogStateVars, vWeights, mShocks)

	function myrewardfunc(vStateVars, vChoiceVars)
		(K, z) = vStateVars
		Kprime = vChoiceVars[1]
		capx = Kprime - (1-δ)*K
		return K^α * exp(z) - capx - γ/2*(capx/K- δ)^2 * K
	end

	function myrewardfunc(Output::Float64, K::Float64, Kprime::Float64)
	    capx = Kprime - (1-δ)K
	    out = (1-τ)*(Output - F*K - γ/2*(capx/K- δ)^2 * K) - (1-κ*(capx<0))*capx + τ * δ * K
	    # dont have condition on F in here!!!!
	end

    function mytransfunc(method::Type{separable}, vExogState, vShocksss)
        # @unpack ρ , σ = p.params
        z = vExogState[1]
        zprime  = ρ*z + σ * vShocksss[1];
        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
        return  inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
    end

	mygrossprofits(vStateVars::Vector{Float64}) = vStateVars[1]^α * exp(vStateVars[2])

    DiscreteDynamicProblem(myrewardfunc, mytransfunc, mygrossprofits,
		params,
		tStateVectors, tChoiceVectors,
        bEndogStateVars, vWeights, mShocks)
end
