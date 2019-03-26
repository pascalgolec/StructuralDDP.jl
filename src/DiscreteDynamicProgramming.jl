module DiscreteDynamicProgramming

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, gridmake, DiscreteDP, VFI
	import QuantEcon: solve # import because will extend
	using Interpolations
	using Parameters
	using NLsolve # to find steady state

	using SparseArrays
	using LinearAlgebra

	export createmodel, solve, drawshocks, simulate
	# export NeoClassicalSimple

	# using DiscreteDynamicModels

	# model definition and parameters
	abstract type DiscreteDynamicModel end
	const DDM = DiscreteDynamicModel

	# abstract type NumChoiceVar end
	# abstract type SingleChoiceVar <: NumChoiceVar end
	# abstract type TwoChoiceVar <: NumChoiceVar end

	# need this for transfunc
	abstract type DDMIntegrationDimension end
	const DDMIntDim = DDMIntegrationDimension
	abstract type StateAction <: DDMIntDim end
	abstract type separable <: DDMIntDim end
	abstract type intermediate <: DDMIntDim end
	const SA = StateAction

	include("helpfunctions.jl")

	# type constructor
	include("DiscreteDynamicProblem.jl")

	# Neoclassical
	include("models/Neoclassical.jl")
	# include("models/Neoclassical/constructor.jl")
	# include("models/Neoclassical/rewardfunc.jl")
	# include("models/Neoclassical/transfunc.jl")
	# include("models/Neoclassical/initialize.jl")

	#Intangible
	include("models/Intangible.jl")
	# include("models/Intangible/constructor.jl")
	# include("models/Intangible/transfunc.jl")
	# include("models/Intangible/rewardfunc.jl")
	# include("models/Intangible/initialize.jl")

	createmodel(sym_model::Symbol; kwargs...) = createmodel(eval(sym_model); kwargs...)

	# Solver
	include("solve/constructors.jl")
	include("solve/main.jl")
	include("solve/rewardmatrix.jl")
	include("solve/transitionmatrix.jl")
	include("solve/singlechoicevar.jl")
	include("solve/twochoicevar.jl")
	include("solve/state_action.jl")
	include("solve/initialendogstatevars.jl")

	# simulator
	include("simulate/shocks.jl")
	include("simulate/simulate.jl")

end
