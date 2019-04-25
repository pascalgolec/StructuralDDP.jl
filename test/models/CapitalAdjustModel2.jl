using NLsolve, Distributions

@doc raw"""
	CapitalAdjustModel2(; kwargs...)

Construct capital investment model with two types of capital and one
productivity shock and convex adjustment costs. The recursive formulation is:

```math
\begin{equation*}
V(K,N,z) = \max_{a_K,a_N} K^\alpha e^z - i_K K - \frac{\gamma_K}{2} i_K^2 K - i_N N - \frac{\gamma_N}{2} i_N^2 N + \beta \mathbb{E} V(K',N',z') \\
i_K \equiv \frac{a_K}{K} - (1-\delta_K) \\
i_N \equiv \frac{a_N}{N} - (1-\delta_N) \\
\text{where } K' = a_K \\
N' = a_N \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
\end{equation*}
```

The firm solves the following problem to initialize at t=0:

```math
\begin{equation*}
V_0(z_0) = \max_{K_0,N_0} V(K_0,N_0,z_0) - (1 + C0) (K_0 + N_0) \\
\text{where } z_0 = \sqrt{\frac{\sigma^2}{1-\rho^2}} ε_0 \\
ε_0 \sim \mathcal{N}(0,1)
\end{equation*}
```

Example:

```julia
using DiscreteDynamicProgramming
prob = CapitalAdjustModel2(nK=30, nN=25, nz=5, ρ=0.5, σ=0.3, δ_K=0.1)
```

""" CapitalAdjustModel2
function CapitalAdjustModel2(;
	β   = 0.90,
	α  = 0.67,
	η  = 0.6,
	δ_K  = 0.1,
	δ_N  = 0.2,
	γ_K  = 1.,
	γ_N  = 2.,
	ρ  = 0.6,
	σ  = 0.3,
	C0  = 0.4,

	nK = 50,
	nN = 50,
	nz = 3,

	# options for testing
	intdim = :Separable_ExogStates,
	initialize_exact = true)


	# initial standard deviation of transitiory shock
    stdz = sqrt(σ^2/(1-ρ^2))
    minz = -3*stdz
    maxz =  3*stdz
    vz = collect(LinRange(minz, maxz, nz))

	function steadystate()
        # check https://www3.nd.edu/~esims1/investment_notes.pdf for details

        function f!(F,x)
            (K, N, I_K, I_N, V_K, V_N) = x
            K = abs(K)
            N = abs(N)

            F[1] = I_K - δ_K * K
            F[2] = I_N - δ_N * N
            F[3] = - 1 - γ_K * I_K/(K+N) + β*V_K
			F[4] = - (1 + γ_N * I_N/(K+N)) + β*V_N
            adjcosts = γ_K/2 * (I_K/(K+N))^2 + γ_N/2 * (I_N/(K+N))^2
            F[5] = - V_K + β*( (exp(stdz^2/2) * α*η*K^(α*η-1)*N^((1-α)*η) +
                                adjcosts) + (1-δ_K)*V_K)
            F[6] = - V_N + β*( (exp(stdz^2/2) * (1-α)*η*K^(α*η)*N^((1-α)*η-1) +
                                adjcosts) + (1-δ_N)*V_N)
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

	function rewardfunc_partial(vStates)
		(K, N, z) = vStates
		return (K^α * N^(1-α))^η * exp(z)
	end

	function adjustcosts(Kprime, K, Nprime, N)

	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N

	    Ktot = K + N

	    adj_K = γ_K/2*(capx/Ktot)^2 * Ktot
	    adj_N = γ_N/2*(intx/Ktot)^2 * Ktot

	    return adj_K + adj_N
	end

	function reward(vStates, vChoices)
		(K, N, z) = vStates
		(Kprime, Nprime) = vChoices

		capx = Kprime - (1-δ_K)*K
		intx = Nprime - (1-δ_N)*N

		return rewardfunc_partial(vStates) - capx - intx -
									adjustcosts(Kprime, K, Nprime, N)
	end

	function reward(Output, vStates, vChoices)
		K, N = vStates
		Kprime, Nprime = vChoices
	    capx = Kprime - (1-δ_K)*K
	    intx = Nprime - (1-δ_N)*N
	    adj = adjustcosts(Kprime, K, Nprime, N)
	    return Output - adj - intx - capx
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

		initializationproblem(value, vChoices) = value - (1+C0)*(vChoices[1] + vChoices[2])

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
		reward,
        transfunc,
		Normal(), # distribution of shocks: standard normal
		β;
		intdim = intdim,
        rewardfunc_partial = rewardfunc_partial,
        initializationproblem = initializationproblem,
        initializefunc = initialize,
		shockdist_initial = shockdist_initial,
		tChoiceVectorsZero = tChoiceVectorsZero,
        )
end
