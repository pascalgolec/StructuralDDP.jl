function NeoClassicalSimple(;
	α = 0.67,
	β = 0.9,
	ρ = 0.6,
	σ = 0.3,
	δ = 0.15,
	γ = 2.0,
	F = 0.01,
	C0 = 0.6,
	π = 0.,
	τ = 0.35,
	κ = 0.,

	nK = 100,
	nz = 20,

	intdim = :Separable_ExogStates)


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

    tStateVectors = (vK, vz) # tuple of basis vectors

	function myrewardfunc(vStateVars, choice)
		(K, z) = vStateVars
		Kprime = choice
		capx = Kprime - (1-δ)*K
		action = (Kprime != K)
		oibdp = K^α * exp(z) - F*K*action - γ/2*(capx/K- δ)^2 * K
		return oibdp*(1-τ) + τ * δ * K - (1-κ*(capx<0))*capx
	end

	# including output, i.e. partial rewardfunc
	function myrewardfunc(Output, vStateVars, Kprime)
        K = vStateVars[1]
	    capx = Kprime - (1-δ)K
	    out = (1-τ)*(Output - F*K - γ/2*(capx/K- δ)^2 * K) - (1-κ*(capx<0))*capx + τ * δ * K
	end

	mygrossprofits(vStateVars) = vStateVars[1]^α * exp(vStateVars[2])

	if intdim == :All

	    transfunc = function mytransfunc(vStates, K, Shock)
	        z = vStates[2]
	        zprime  = ρ*z + σ * Shock;
	        return K, zprime
	    end
		tChoiceVectors = (vK,)

	elseif intdim == :Separable

	    transfunc = function mytransfunc1(vStates, vChoices, Shock)
	        z = vStates[2]
	        zprime  = ρ*z + σ * Shock;
	        return zprime
	    end
		tChoiceVectors = (1,)

	elseif intdim == :Separable_States
	    transfunc = function mytransfunc2(vStates, Shock)
	        z = vStates[2]
	        zprime  = ρ*z + σ * Shock;
	        return zprime
	    end
		tChoiceVectors = (1,)

	elseif intdim == :Separable_ExogStates
	    transfunc = function mytransfunc3(ExogState, shock)
	        z = ExogState
	        zprime  = ρ*z + σ * shock;
	        return zprime
	    end
		tChoiceVectors = (1,)
	else

		error("$intdim wrong intdim")
	end


	tChoiceVectorsZero = (1,)

	initializationproblem(value, choice) = value - (1 + (1-β)/β + C0) * choice

	function initialize(shock)
	    z0 = shock * sqrt(σ^2 / (1-ρ^2))
		return z0
	end

	DDP(		tStateVectors,
	            tChoiceVectors,
	            myrewardfunc,
	            transfunc,
				Normal(), # give distribution of shocks: standard normal
				β;
				intdim = intdim,
	            rewardfunc_partial = mygrossprofits,
	            initializationproblem = initializationproblem,
	            initializefunc = initialize,
				tChoiceVectorsZero = tChoiceVectorsZero,
	            )

	# DDP(
	#
	# 				tStateVectors,
	# 				tChoiceVectors,
	# 				myrewardfunc,
	# 				transfunc,
	# 				Separable_ExogStates,
	# 				Normal(),
	# 				β,
	# 				nothing)
end
