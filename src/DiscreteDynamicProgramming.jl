module DiscreteDynamicProgramming

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, gridmake
	using Interpolations
	using Parameters
	using NLsolve # to find steady state

	using SparseArrays
	using LinearAlgebra

	export createmodel, solve
	export NeoClassicalSimple

	# model definition and parameters
	abstract type DiscreteDynamicModel end
	const DDM = DiscreteDynamicModel
	abstract type SingleChoiceVar <: DDM end

	abstract type DDMIntegrationDimension end
	const DDMIntDim = DDMIntegrationDimension
	abstract type StateAction <: DDMIntDim end
	abstract type separable <: DDMIntDim end
	abstract type intermediate <: DDMIntDim end
	const SA = StateAction

	include("helpfunctions.jl")

	# Neoclassical
	include("models/Neoclassical/constructor.jl")
	include("models/Neoclassical/initializationproblem.jl")
	include("models/Neoclassical/initialize.jl")
	include("models/Neoclassical/outputfunc.jl")
	include("models/Neoclassical/rewardfunc.jl")
	include("models/Neoclassical/transfunc.jl")

	# Solver
	include("solve/constructors.jl")
	include("solve/main.jl")
	include("solve/rewardmatrix.jl")
	include("solve/transitionmatrix.jl")
	include("solve/singlechoicevar.jl")
	include("solve/initialendogstatevars.jl")

end
