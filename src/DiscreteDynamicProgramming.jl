module DiscreteDynamicProgramming

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, qnwlogn, gridmake, DiscreteDP, VFI
	using Interpolations, Distributions, TreeViews, DataFrames
	using SparseArrays, LinearAlgebra
	using Test

	import QuantEcon: solve
	import DataFrames: DataFrame
	import Base: Array

	using DocStringExtensions

	export DiscreteDynamicProblem, DDP
	export createmodel
	export solve, transitionmatrix, value, policy, compare
	export drawshocks, simulate, policy, states, value, DataFrame

	include("utils.jl")

	# type constructor
	include("problem/Transition.jl")
	include("problem/DiscreteDynamicProblem.jl")
	include("problem/DiscreteDynamicProblem_interface.jl")

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
	include("simulate/constructor.jl")
	include("simulate/shocks.jl")
	include("simulate/simulate.jl")

	# Models for testing
	include("../test/models/Neoclassical_user.jl")
	include("../test/models/Intangible_user.jl")
	createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)

end
