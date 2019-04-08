module DiscreteDynamicProgramming

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, qnwlogn, gridmake, DiscreteDP, VFI
	import QuantEcon: solve # import because will extend
	using Interpolations
	using Parameters
	using NLsolve # to find steady state
	using TreeViews
	using Distributions

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
	abstract type All <: DDMIntDim end
	abstract type Separable <: DDMIntDim end
	abstract type Separable_States <: DDMIntDim end
	abstract type Separable_ExogStates <: DDMIntDim end
	const Separable_Union = Union{Separable, Separable_States, Separable_ExogStates}
	# const SA = StateAction

	include("utils.jl")

		# type constructor
	include("DiscreteDynamicProblem.jl")
	# include("DiscreteDynamicProblem_interface.jl")

	# Models
	include("models/Neoclassical_user.jl")
	# include("models/Intangible.jl")
	# include("models/Intangible.jl")

	# for using createmodel syntax:
	createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)
	# e.g. createmodel(:NeoClassicalSimple; nK = 30, nz=15)

	# # Solver
	include("solve/constructors.jl")
	include("solve/main.jl")
	include("solve/rewardmatrix.jl")
	include("solve/transitionmatrix.jl")
	include("solve/singlechoicevar.jl")
	# include("solve/twochoicevar.jl")
	include("solve/state_action.jl")
	# include("solve/initialendogstatevars.jl")

	# # simulator
	# include("simulate/shocks.jl")
	# include("simulate/simulate.jl")

end
