# Model: Peters and Taylor 2017
# profits = (Kphy^α Kint^(1-α))^γ
# law of motion: 	Kphy' = (1-δphy)Kphy + Iphy
# 					Kint' = (1-δint)Kint + Iint
# adjustment costs:
#   - only convex, not fixed
@with_kw struct IntangibleParams{R<:Real, I<:Int} <: ModelParams
	β ::R  = 0.96
	α ::R = 0.67
	γ ::R = 0.6
	δ_K ::R = 0.1
	δ_N ::R = 0.2
	a_K ::R = 1.
	a_N ::R = 2.
	ρ ::R = 0.6
	σ ::R = 0.3
	C0 ::R = 0.4
	τ ::R = 0.3

	nK::I = 50
	nN::I = 50
	nz::I = 3

	nShocks::I = 3
	nPeriods::I = 40
	nFirms::I = 100

end

function Intangible(; kwargs...)

	params = IntangibleParams(;kwargs...)
    @unpack_IntangibleParams params

	# initial standard deviation of transitiory shock
    stdz = sqrt(σ^2/(1-ρ^2))
    minz = -3*stdz
    maxz =  3*stdz
    vz = collect(LinRange(minz, maxz, nz))

	function steadystate()
        # check https://www3.nd.edu/~esims1/investment_notes.pdf

        function f!(F,x)
            (K, N, I_K, I_N, V_K, V_N) = x
            K = abs(K)
            N = abs(N)

            F[1] = I_K - δ_K * K
            F[2] = I_N - δ_N * N
            F[3] = - 1 - a_K * I_K/(K+N)*(1-τ) + β*V_K
			F[4] = - (1 + a_N * I_N/(K+N))*(1-τ) + β*V_N
            adjcosts = a_K/2 * (I_K/(K+N))^2 + a_N/2 * (I_N/(K+N))^2
            F[5] = - V_K + β*( (exp(stdz^2/2) * α*γ*K^(α*γ-1)*N^((1-α)*γ) +
                                adjcosts)*(1-τ) + τ*δ_K + (1-δ_K)*V_K)
            F[6] = - V_N + β*( (exp(stdz^2/2) * (1-α)*γ*K^(α*γ)*N^((1-α)*γ-1) +
                                adjcosts)*(1-τ) + (1-δ_N)*V_N)
        end

        out = nlsolve(f!, [ 10.; 10.; 0.1; 0.1; 1.1; 1.1])

        # println("K = ", abs(out.zero[1]))
        # println("logK = ", log(abs(out.zero[1])))
        # println("N = ", abs(out.zero[2]))
        # println("logN = ", log(abs(out.zero[2])))
        # println("I_K = ", out.zero[3])
        # println("I_N = ", out.zero[4])

        return (abs(out.zero[1]), abs(out.zero[2]))

    end

    (K_ss, N_ss) = steadystate()

    # K_ss = (α*γ * exp(stdz^2/2)/ ((1-β)/β + δ_K)) ^ (1/(1-α*γ))
    # minK = exp(0.7 * log(K_ss))
    # maxK = exp(1.25 * log(K_ss))
	minK = 0.1 * K_ss
    maxK = 3.5 * K_ss
    vK   = exp.(collect(LinRange(log(minK), log(maxK), nK)))

    # N_ss = 123456.
    # minN = 0.5
    # maxN = 15.
    # N_ss = ((1-α)*γ * exp(stdz^2/2)/ ((1-β)/β + δ_N)) ^ (1/(1-(1-α)*γ))
    # minN = exp(0.7 * log(N_ss))
    # maxN = exp(1.25 * log(N_ss))
	minN = 0.1 * N_ss
    maxN = 3.5 * N_ss
 	vN   = exp.(collect(LinRange(log(minN), log(maxN), nN)))

    # vMin = [minK, minN, minz]
    # vMax = [maxK, maxN, maxz]

    tStateVectors = (vK, vN, vz)
    tChoiceVectors = (vK, vN)

    bEndogStateVars = [true, true, false]

    # shocks
    vShocks, vWeights = qnwnorm(nShocks,0,1)
    mShocks = Array(vShocks') # need transpose because want row vector

	# Intangible(params, tStateVectors, tChoiceVectors,
    #     bEndogStateVars, vWeights, mShocks)

	function mygrossprofits(vStateVars::Vector{Float64})
		(K, N, z) = vStateVars
		return (K^α * N^(1-α))^γ * exp(z)
	end

	function myadjustcosts(Kprime::Float64,
	                     K::Float64,
						 Nprime::Float64,
						 N::Float64)

	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N

	    Ktot = K + N

	    adj_K = a_K/2*(capx/Ktot)^2 * Ktot
	    adj_N = a_N/2*(intx/Ktot)^2 * Ktot

	    return adj_K + adj_N
	end

	function myrewardfunc(vStateVars::Vector{Float64}, vChoices::Vector{Float64})
		(K, N, z) = vStateVars
		(Kprime, Nprime) = vChoices

		capx = Kprime - (1-δ_K)*K
		intx = Nprime - (1-δ_N)*N

		oibdp = mygrossprofits(vStateVars) -
									myadjustcosts(Kprime, K, Nprime, N) - intx
		return oibdp*(1-τ) + τ*δ_K*K - capx
	end

	function myrewardfunc(Output::Float64, K::Float64, Kprime::Float64,
	                                     N::Float64, Nprime::Float64)
	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N

	    adj = myadjustcosts(Kprime, K, Nprime, N)

	    return (Output - adj - intx)*(1-τ) - capx + δ_K*τ*K
	end


    function mytransfunc(method::Type{separable}, vExogState::Vector{Float64}, shock::Float64)
        # @unpack ρ , σ = p.params
        z = vExogState[1]
        zprime  = ρ*z + σ * shock;
        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
        return  inbounds(zprime, tStateVectors[3][1], tStateVectors[3][end])
    end
    function mytransfunc(method::Type{separable}, vExogState::Vector{Float64}, vShocksss::Vector{Float64})
        # @unpack ρ , σ = p.params
        z = vExogState[1]
        zprime  = ρ*z + σ * vShocksss[1];
        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
        return  inbounds(zprime, tStateVectors[3][1], tStateVectors[3][end])
    end

	# this is superflous, only for testing
	function mytransfunc(method::Type{intermediate}, vState::Vector{Float64}, shock::Float64)
        z = vState[3]
        zprime  = ρ*z + σ * shock
        return  inbounds(zprime, tStateVectors[3][1], tStateVectors[3][end])
    end
	function mytransfunc(method::Type{intermediate}, vState::Vector{Float64}, vShocksss::Vector{Float64})
        z = vState[3]
        zprime  = ρ*z + σ * vShocksss[1]
        return  inbounds(zprime, tStateVectors[3][1], tStateVectors[3][end])
    end

	function mytransfunc(method::Type{SA}, vState::Vector{Float64}, vChoice::Vector{Float64}, vShocksss::Vector{Float64})
        # @unpack ρ , σ = p.params
		# @code_warntype vState[2]
        z = vState[3]
		# @code_warntype ρ*z + σ * vShocksss[1]
        zprime  = ρ*z + σ * vShocksss[1]
        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
        return  [vChoice[1], vChoice[2], inbounds(zprime, tStateVectors[3][1], tStateVectors[3][end])]
    end

	function initializationproblem(value::Float64, choiceone::Float64, choicetwo::Float64)
	    K = choiceone
	    N = choicetwo
	    # want to choose K to maximize V[K,N,z] - (1+p.C0)*K - (1+p.C0)*N
	    return value - (1+C0)*K - (1+C0)*N
	end

	function initialize(vShocks::Array{Float64, 1}, itp_K0, itp_N0)

	    z0 = vShocks[1] * sqrt(σ^2 / (1-ρ^2))
		z0 = inbounds(z0, tStateVectors[3][1], tStateVectors[3][end])

		K0 = itp_K0(z0)
		K0 = inbounds(K0, tStateVectors[1][1], tStateVectors[1][end])

		N0 = itp_N0(z0)
		N0 = inbounds(N0, tStateVectors[2][1], tStateVectors[2][end])

		return [K0, N0, z0]
	end


	DiscreteDynamicProblem(	params,
	            myrewardfunc,
	            mytransfunc,
	            separable,
	            tStateVectors,
	            tChoiceVectors,
	            vWeights,
	            mShocks;
	            bEndogStateVars = bEndogStateVars,
	            grossprofits = mygrossprofits,
	            initializationproblem = initializationproblem,
	            initializefunc = initialize,
	            )
end
