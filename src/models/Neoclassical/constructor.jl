
####################
# NeoClassical
####################

@with_kw struct NeoClassicalSimpleParams{R<:Real, I<:Int64, B<:Bool}
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

    rewardmat::Symbol = :nobuild
    intdim::Symbol = :separable
    monotonicity::B = true
    concavity::B = true
end

struct NeoClassicalSimple <: SingleChoiceVar
    # parameters that user supplies
    params::NeoClassicalSimpleParams{Float64, Int64, Bool} # important to specify here for type stability

    # K_ss::Float64

    # parameters that are a function of inputs or always the same depending on model
    # vMin::Vector{Float64}
    # vMax::Vector{Float64}

    tStateVectors::NTuple{2, Vector{Float64}}
    tChoiceVectors::NTuple{1, Vector{Float64}}

    # which state variables are endogenous
    bEndogStateVars::Vector{Bool}

    # quadrature for calculating expectations
    vWeights::Vector{Float64}
    mShocks::Array{Float64,2}
end

createmodel(sym_model; kwargs...) = createmodel(eval(sym_model); kwargs...)

function createmodel(model::Type{NeoClassicalSimple}; kwargs...)

    # all userprovided parameters go in here
    params = NeoClassicalSimpleParams(;kwargs...)

    # need to unpack to get default values if not provided
    @unpack_NeoClassicalSimpleParams params

    nNodes = collect((nK, nz))

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

    # NeoClassicalSimple(params, K_ss, vMin, vMax, tStateVectors, tChoiceVectors,
    #     bEndogStateVars, vWeights, mShocks)
    NeoClassicalSimple(params, tStateVectors, tChoiceVectors,
        bEndogStateVars, vWeights, mShocks)
end
