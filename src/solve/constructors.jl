
# more types for algorhythm
# abstract type capitalaction end
# abstract type active <: capitalaction end
# abstract type passive <: capitalaction end
# # only use capital action for multiplechoicevar. will try to code without

"""ItpArray allows the array to be interpolated. Remove when Interpolations.jl
supports array indexing with `[]`."""
struct ItpArray{T,N,Tarr<:Array{T,N},
		Titp<:AbstractInterpolation{T,N}} <: AbstractArray{T,N}
	arr::Tarr
	itp::Titp
end
# overload array functions
Base.size(A::ItpArray) = size(A.arr)
Base.getindex(A::ItpArray, args...) = getindex(A.arr,args...)
Base.setindex!(A::ItpArray, args...) = setindex!(A.arr, args...)
# overload interpolator
function (A::ItpArray)(x::Vararg{Number,N}) where {T,N}
	return A.itp(x...) #A.itp()(x) # size(A.itp)
end
function itparray(A::AbstractArray{T,N}, knots::NTuple{N,AbstractVector{T}}) where {T,N}
	itp = interpolate(knots, A, Gridded(Linear()))
	itp = extrapolate(itp, Flat())
	ItpArray{T,N,typeof(A),typeof(itp)}(A,itp)
end

abstract type AbstractDDPSolution{NS,NC} end

"""Solution object.
Parameters:
    - NS = number of states variables
    - NC = number of choice variables"""
struct DDPSolution{NS,NC,valueType,policyType} <: AbstractDDPSolution{NS,NC}

    prob::DDP
    value::valueType
    policy::NTuple{NC,policyType}

end

"""Includes solution for exact initialization.
Parameters:
    - NE0 = number of exogenous states variables at t=0
    - NC0 = number of choice variables at t=0"""
struct DDPSolutionZero{NS,NC,NE0,NC0,valueType,policyType,
		value0Type,policy0Type} <: AbstractDDPSolution{NS,NC}

	prob::DDP
	value::valueType
	policy::NTuple{NC,policyType}
	value0::value0Type
	policy0::NTuple{NC0,policy0Type}

end

function createsolution(p::DDP, meshValFun::Array{T,NS},
                                tmeshPolFun::NTuple{NC,Array{T,NS}}) where {T,NS,NC}

	value = itparray(meshValFun, p.tStateVectors)
  	policy = [itparray(pol, p.tStateVectors) for pol in tmeshPolFun]
  	policy = tuple(policy...)

   if p.options.initialize == nothing
      return DDPSolution{NS,NC,typeof(value),typeof(policy[1])}(p, value, policy)
   else
      meshValFunZero, tmeshPolFunZero = initialendogstatevars(p, meshValFun)
		tExogStateVectorsZero = getnonchoicevarszero(p)
		NE0 = length(tExogStateVectorsZero)
		NC0 = length(p.options.initialize.tChoiceVectorsZero)
		value0 = itparray(meshValFunZero, tExogStateVectorsZero)
		policy0 = [itparray(pol, tExogStateVectorsZero) for pol in tmeshPolFunZero]
		policy0 = tuple(policy0...)
      return DDPSolutionZero{NS,NC,NE0,NC0,typeof(value),typeof(policy[1]),
			typeof(value0),typeof(policy0[1])}(p, value, policy, value0, policy0)
   end
end

value(sol::AbstractDDPSolution) = sol.value
policy(sol::AbstractDDPSolution) = sol.policy
policy(sol::AbstractDDPSolution{NS,1}) where NS = sol.policy[1]

value0(sol::DDPSolutionZero) = sol.value0
policy0(sol::DDPSolutionZero) = sol.policy0
policy0(sol::DDPSolutionZero{NS,1}) where NS  = sol.policy0[1]
