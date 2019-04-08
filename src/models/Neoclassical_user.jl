function NeoClassicalSimple(; 	α = 0.67,
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

								nPeriods = 120,
								nFirms = 1000,

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

    # is incorrect
    # K_ss_analytical(z) = (α * (1-τ) * exp(z)/ ((1-β)/β*(1+γ*δ) + δ)) ^ (1/(1-α))
    # @show log(K_ss_analytical(minz))
    # @show log(K_ss_analytical(maxz))

    # vMin = [exp(minK_log), minz]
    # vMax = [exp(maxK_log), maxz]

    tStateVectors = (vK, vz) # tuple of basis vectors

	# tChoiceVectors = (vK,) # important to add the comma, otherwise not a tuple
    # tChoiceVectors = (1,)

    # bEndogStateVars = [true, false]

	# somehow this below is not type stable!
	# myrewardfunc(vStateVars::Vector{Float64}, vChoices::Vector{Float64}) =
	#     myrewardfunc(vStateVars, vChoices[1])

	function myrewardfunc(vStateVars, vChoices)
		(K, z) = vStateVars
		# Kprime = vChoiceVars[1]
		Kprime = vChoices[1]
		capx = Kprime - (1-δ)*K
		action = (Kprime != K)
		oibdp = K^α * exp(z) - F*K*action - γ/2*(capx/K- δ)^2 * K
		return oibdp*(1-τ) + τ * δ * K - (1-κ*(capx<0))*capx
	end

	# including output, i.e. partial rewardfunc
	# function myrewardfunc(Output::Float64, K::Float64, Kprime::Float64)
	function myrewardfunc(Output, K, Kprime)
	    capx = Kprime - (1-δ)K
	    out = (1-τ)*(Output - F*K - γ/2*(capx/K- δ)^2 * K) - (1-κ*(capx<0))*capx + τ * δ * K
	    # dont have condition on F in here!!!!
	end

	mygrossprofits(vStateVars) = vStateVars[1]^α * exp(vStateVars[2])

	if intdim == :All

	    transfunc = function mytransfunc(vStates, vChoices, vShocksss)
	        # @unpack ρ , σ = p.params
	        z = vStates[2]
	        zprime  = ρ*z + σ * vShocksss[1];
	        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
	        return vChoices[1], inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
	    end
		tChoiceVectors = (vK,)

	elseif intdim == :Separable

	    transfunc = function mytransfunc1(vStates, vChoices, vShocksss)
	        # @unpack ρ , σ = p.params
	        z = vStates[2]
	        zprime  = ρ*z + σ * vShocksss[1];
	        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
	        return inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
	    end
		tChoiceVectors = (1,)

	elseif intdim == :Separable_States
	    transfunc = function mytransfunc2(vStates, vShocksss)
	        # @unpack ρ , σ = p.params
	        z = vStates[2]
	        zprime  = ρ*z + σ * vShocksss[1];
	        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
	        return inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
	    end
		tChoiceVectors = (1,)

	elseif intdim == :Separable_ExogStates
	    transfunc = function mytransfunc3(vExogState, vShocksss)
	        # @unpack ρ , σ = p.params
	        z = vExogState[1]
	        zprime  = ρ*z + σ * vShocksss[1];
	        # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
	        return inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
	    end
		tChoiceVectors = (1,)
	else
		error("$intdim wrong intdim")
	end

	# # specification is with endogchoicevars but SA
	# function mytransfunc(vState::Vector{Float64}, vChoice::Vector{Float64}, vShocksss::Vector{Float64})
    #     # @unpack ρ , σ = p.params
	# 	# @code_warntype vState[2]
    #     z = vState[2]
	# 	# @code_warntype ρ*z + σ * vShocksss[1]
    #     zprime  = ρ*z + σ * vShocksss[1]
    #     # return  inbounds(zprime, p.tStateVectors[2][1], p.tStateVectors[2][end])
    #     return  inbounds(zprime, tStateVectors[2][1], tStateVectors[2][end])
    # end

	initializationproblem(value::Float64, K::Float64) =
		value - (1 + (1-β)/β + C0) * K

	function initialize(dShock::AbstractArray{Float64, 1}, itp_K0)
	    z0 = dShock[1] * sqrt(σ^2 / (1-ρ^2))
		z0 = inbounds(z0, tStateVectors[2][1], tStateVectors[2][end])
	    K0 = itp_K0(z0)
		K0 = inbounds(K0, tStateVectors[1][1], tStateVectors[1][end])
		return [K0, z0]
	end

	createDiscreteDynamicProblem(
				β,
	            myrewardfunc,
	            transfunc,
	            tStateVectors,
	            tChoiceVectors,
				Normal(); # give distribution of shocks: standard normal
				intdim = intdim,
	            # bEndogStateVars = bEndogStateVars,
	            grossprofits = mygrossprofits,
	            initializationproblem = initializationproblem,
	            initializefunc = initialize,
	            )
end
