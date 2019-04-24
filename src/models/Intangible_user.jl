function Intangible(;
	β   = 0.90,
	α  = 0.67,
	γ  = 0.6,
	δ_K  = 0.1,
	δ_N  = 0.2,
	a_K  = 1.,
	a_N  = 2.,
	ρ  = 0.6,
	σ  = 0.3,
	C0  = 0.4,
	τ  = 0.3,

	nK = 50,
	nN = 50,
	nz = 3,

	intdim = :Separable_ExogStates,
	initialize_exact = true)


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
        return (abs(out.zero[1]), abs(out.zero[2]))
    end

    (K_ss, N_ss) = steadystate()

	minK = 0.1 * K_ss
    maxK = 3.5 * K_ss
    vK   = exp.(collect(LinRange(log(minK), log(maxK), nK)))

	minN = 0.1 * N_ss
    maxN = 3.5 * N_ss
 	vN   = exp.(collect(LinRange(log(minN), log(maxN), nN)))

    tStateVectors = (vK, vN, vz)

	function mygrossprofits(vStateVars)
		(K, N, z) = vStateVars
		return (K^α * N^(1-α))^γ * exp(z)
	end

	function myadjustcosts(Kprime, K, Nprime, N)

	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N

	    Ktot = K + N

	    adj_K = a_K/2*(capx/Ktot)^2 * Ktot
	    adj_N = a_N/2*(intx/Ktot)^2 * Ktot

	    return adj_K + adj_N
	end

	function myrewardfunc(vStateVars, vChoices)
		(K, N, z) = vStateVars
		(Kprime, Nprime) = vChoices

		capx = Kprime - (1-δ_K)*K
		intx = Nprime - (1-δ_N)*N

		oibdp = mygrossprofits(vStateVars) -
									myadjustcosts(Kprime, K, Nprime, N) - intx
		return oibdp*(1-τ) + τ*δ_K*K - capx
	end

	function myrewardfunc(Output, vStateVars, vChoices)
		K, N = vStateVars
		Kprime, Nprime = vChoices
	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N
	    adj = myadjustcosts(Kprime, K, Nprime, N)
	    return (Output - adj - intx)*(1-τ) - capx + δ_K*τ*K
	end

	if intdim == :All
	    transfunc = function mytransfunc(vStates, vChoices, Shock)
	        z = vStates[3]
	        zprime  = ρ*z + σ * Shock
	        return vChoices[1], vChoices[2], zprime
	    end
		tChoiceVectors = (vK, vN)

	elseif intdim == :Separable
	    transfunc = function mytransfunc1(vStates, vChoices, Shock)
	        z = vStates[3]
	        zprime  = ρ*z + σ * Shock
	        return zprime
	    end
		tChoiceVectors = (1, 2)

	elseif intdim == :Separable_States
	    transfunc = function mytransfunc2(vStates, Shock)
	        z = vStates[3]
	        zprime  = ρ*z + σ * Shock
	        return zprime
	    end
		tChoiceVectors = (1, 2)

	elseif intdim == :Separable_ExogStates
	    transfunc = function mytransfunc3(z, Shock)
	        zprime  = ρ*z + σ * Shock
	        return zprime
	    end
		tChoiceVectors = (1, 2)

	else
		error("$intdim wrong intdim")
	end

	if initialize_exact
		tChoiceVectorsZero = (1,2)

		initializationproblem(value, vChoices) = value - (1+C0)*vChoices[1] - (1+C0)*vChoices[2]

		function initialize(shock)
		    z0 = shock * sqrt(σ^2 / (1-ρ^2))
			# z0 = inbounds(z0, tStateVectors[3][1], tStateVectors[3][end])
			return z0
		end
		shockdist_initial = Normal()
	else
		tChoiceVectorsZero = nothing
		initializationproblem = nothing
		initialize = nothing
		shockdist_initial = nothing
	end

	DDP(tStateVectors,
		tChoiceVectors,
		myrewardfunc,
        transfunc,
		Normal(), # distribution of shocks: standard normal
		β;
		intdim = intdim,
        rewardfunc_partial = mygrossprofits,
        initializationproblem = initializationproblem,
        initializefunc = initialize,
		shockdist_initial = shockdist_initial,
		tChoiceVectorsZero = tChoiceVectorsZero,
        )
end
