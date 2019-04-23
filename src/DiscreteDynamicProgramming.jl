module DiscreteDynamicProgramming

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, qnwlogn, gridmake, DiscreteDP, VFI
	using Interpolations, Distributions, TreeViews, DataFrames
	using SparseArrays, LinearAlgebra
	using NLsolve # to find steady state
	using Test

	import QuantEcon: solve
	import DataFrames: DataFrame
	import Base: Array

	using DocStringExtensions

	export DiscreteDynamicProblem
	export createmodel
	export solve, transitionmatrix, value, policy
	export drawshocks, simulate, policy, states, value
	export compare

	# using DiscreteDynamicModels

	include("utils.jl")

	# type constructor
	include("DiscreteDynamicProblem.jl")
	# include("DiscreteDynamicProblem_interface.jl")

	# Models
	include("models/Neoclassical_user.jl")
	include("models/Intangible_user.jl")

	# for using createmodel syntax:
	createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)
	# e.g. createmodel(:NeoClassicalSimple; nK = 30, nz=15)

	# Solver
	include("solve/constructors.jl")
	include("solve/main.jl")
	include("solve/rewardmatrix.jl")
	include("solve/transitionmatrix.jl")
	include("solve/singlechoicevar.jl")
	include("solve/twochoicevar.jl")
	include("solve/state_action.jl")
	include("solve/initialendogstatevars.jl")
	include("solve/compare.jl")

	# simulator
	include("simulate/constructor.jl")
	include("simulate/shocks.jl")
	include("simulate/simulate.jl")

end
