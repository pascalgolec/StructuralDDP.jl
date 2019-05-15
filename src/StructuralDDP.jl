module StructuralDDP

	using BasisMatrices: Basis, SplineParams, nodes, BasisMatrix, Expanded
	using QuantEcon: qnwnorm, qnwunif, qnwlogn, gridmake, DiscreteDP, VFI
	using Interpolations, Distributions, TreeViews, DataFrames
	using SparseArrays, LinearAlgebra
	using DocStringExtensions
	using Parameters

	import QuantEcon: solve
	import DataFrames: DataFrame
	import Base: Array, isapprox

	export DiscreteDynamicProblem, DDP
	export createmodel
	export solve, transitionmatrix, value, policy
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

	# Models for developing: requires NLsolve
	include("../test/models/CapitalAdjustModel.jl")
	include("../test/models/CapitalAdjustModel2.jl")
	createmodel(model::Symbol; kwargs...) = eval(model)(; kwargs...)
	export CapitalAdjustModel, CapitalAdjustModel2, createmodel

end
