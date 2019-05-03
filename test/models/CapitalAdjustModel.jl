using NLsolve, Distributions

@doc raw"""
	CapitalAdjustModel(; kwargs...)

Construct capital investment model with one type of capital and one
productivity shock and convex adjustment costs. The recursive formulation is:

```math
\begin{equation*}
V(K,z) = \max_a K^\alpha e^z - i K - \frac{\gamma}{2} i^2 K + \beta \mathbb{E} V(K', z') \\
i \equiv \frac{a}{K} - (1-\delta) \\
\text{where } K' = a \\
z' = \rho z + \sigma \varepsilon, \quad \varepsilon \sim \mathcal{N}(0,1)
\end{equation*}
```

The firm solves the following problem to initialize at t=0:

```math
\begin{equation*}
V_0(z_0) = \max_{K_0} V(K_0,z_0) - (1 + C0) K_0 \\
\text{where } z_0 = \sqrt{\frac{\sigma^2}{1-\rho^2}} ε_0 \\
ε_0 \sim \mathcal{N}(0,1)
\end{equation*}
```

Example:

```julia
using StructuralDDP
prob = CapitalAdjustModel(nK=150, nz=15, ρ=0.5, σ=0.3, γ=2.)
```

"""
function CapitalAdjustModel(;
	α = 0.67,
	β = 0.9,
	ρ = 0.6,
	σ = 0.3,
	δ = 0.15,
	γ = 2.0,
	C0 = 0.6,
	F = 0.,

	nK = 100,
	nz = 20,

	# options for testing
	intdim = :Separable_ExogStates,
	initialize_exact = true)


    # calculate Amin, Amax
    stdz = sqrt(σ^2/(1-ρ^2))
    minz = -3*stdz
    maxz =  3*stdz
    vz = collect(LinRange(minz, maxz, nz))

    function steadystate(z)
        # check https://www3.nd.edu/~esims1/investment_notes.pdf

        function f2!(F,x)
            (V_K, K_log) = x

            F[1] = - 1 - γ + V_K
            F[2] = V_K - β*( exp(z)*α*exp(K_log)^(α-1) + (1-δ)*V_K)
        end

        return nlsolve(f2!, [1.1; 1.])
    end

    K_ss = exp(steadystate(stdz^2/2).zero[2])
    minK_log = steadystate(-2*stdz).zero[2]
    maxK_log = steadystate(3*stdz).zero[2]
    vK   = exp.(collect(LinRange(minK_log, maxK_log, nK)))

    tStateVectors = (vK, vz) # tuple of basis vectors

	function reward(vStates, vChoices)
	    K, z = vStates
	    Kprime = vChoices[1]
		i = Kprime/K - (1-δ)
		action = Kprime != K
	    return K^α * exp(z) - i*K - F*K*action - γ/2 * i^2 * K
	end

	# including output, i.e. for partial rewardfunc
	function reward(partial_reward, vStates, vChoices)
	    K = vStates[1]
		Kprime = vChoices[1]
	    i = Kprime/K - (1-δ)
		action = Kprime != K
	    return partial_reward - i*K - F*K*action - γ/2 * i^2 * K
	end

	reward_partial(vStates) = vStates[1]^α * exp(vStates[2])

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
	    transfunc = function mytransfunc3(vExogStates, shock)
	        z = vExogStates[1]
	        zprime  = ρ*z + σ * shock;
	        return zprime
	    end
		tChoiceVectors = (1,)
	else

		error("$intdim wrong intdim")
	end


	if initialize_exact
		tChoiceVectorsZero = (1,)

		initializationproblem(value, choice) = value - (1 + (1-β)/β + C0) * choice

		function initialize(shock)
		    z0 = shock * sqrt(σ^2 / (1-ρ^2))
			return z0
		end
		shockdist_initial = Normal()
	else
		tChoiceVectorsZero = nothing
		initializationproblem = nothing
        initialize = nothing
		shockdist_initial = nothing
	end

	if F >0
		function checkwhich(vStatesIndex)
			vChoice = vStatesIndex[1]
			return vChoice
		end
	else
		checkwhich = nothing
	end

	DDP(tStateVectors,
        tChoiceVectors,
        reward,
        transfunc,
		Normal(), # give distribution of shocks: standard normal
		β;
		intdim = intdim,
        rewardfunc_partial = reward_partial,
        initializationproblem = initializationproblem,
        initializefunc = initialize,
		shockdist_initial = shockdist_initial,
		tChoiceVectorsZero = tChoiceVectorsZero,
		get_additional_index = checkwhich,
        )
end
