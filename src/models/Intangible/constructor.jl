# Model:
# profits = (Kphy^α Kint^(1-α))^γ
# law of motion: 	Kphy' = (1-δphy)Kphy + Iphy
# 					Kint' = (1-δint)Kint + Iint
# adjustment costs:
#   - only convex, not fixed
@with_kw struct IntangibleParams{R<:Real, I<:Int64}
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

	nShocks::I = 2
	nPeriods::I = 40
	nFirms::I = 100

	# rewardmat::Symbol = :nobuild
	# intdim::Symbol = :separable
	# monotonicity::Vector{B} = [true, true]
	# concavity::Vector{B} = [true, true]
end

struct Intangible <: TwoChoiceVar
    # parameters that user supplies
    params::IntangibleParams{Float64, Int64}

    tStateVectors::NTuple{3, Vector{Float64}}
    tChoiceVectors::NTuple{2, Vector{Float64}}

    # which state variables are endogenous
    bEndogStateVars::Vector{Bool}

    # quadrature for calculating expectations
    vWeights::Vector{Float64}
    mShocks::Array{Float64,2}
end

function createmodel(model::Type{Intangible}; kwargs...)
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
    mShocks = vShocks' # need transpose because want row vector

	Intangible(params, tStateVectors, tChoiceVectors,
        bEndogStateVars, vWeights, mShocks)
end
