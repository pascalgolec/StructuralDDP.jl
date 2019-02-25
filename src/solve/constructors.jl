
# more types for algorhythm
# abstract type capitalaction end
# abstract type active <: capitalaction end
# abstract type passive <: capitalaction end
# # only use capital action for multiplechoicevar. will try to code without

struct DDPSolution

       meshValFun::Array{Float64}
       tmeshPolFun::NTuple{N,Array{Float64}} where N
       # mV0::Array{Float64}
       # mPolFun0::Array{Float64}

end


function createsolution(p::DDM, meshValFun::Array{Float64},
                                tmeshPolFun::NTuple{N,Array{Float64}}) where
                                N # implies singlechoicevar

        # mV0, mPolFun0 = initialendogstatevars(p, meshValFun)

        # SingleChoiceVarSolution(meshValFun, meshPolFun, mV0, mPolFun0)
        DDPSolution(meshValFun, tmeshPolFun)
end

# abstract type DDPSolution end

# constructors

# struct SingleChoiceVarSolution <: DDPSolution
#
#        meshValFun::Array{Float64}
#        meshPolFun::Array{Float64}
#        # mV0::Array{Float64}
#        # mPolFun0::Array{Float64}
#
# end


# function createsolution(p::SingleChoiceVar, meshValFun::Array{Float64},
#                                             meshPolFun::Array{Float64})
#
#         # mV0, mPolFun0 = initialendogstatevars(p, meshValFun)
#
#         # SingleChoiceVarSolution(meshValFun, meshPolFun, mV0, mPolFun0)
#         SingleChoiceVarSolution(meshValFun, meshPolFun)
# end

# struct TwoChoiceVarolution <: DDPSolution
#
#        meshValFun::Array{Float64}
#        meshPolFun1::Array{Float64}
#        meshPolFun2::Array{Float64}
#        mV0::Array{Float64}
#        mPolFunOne0::Array{Float64}
#        mPolFunTwo0::Array{Float64}
#
# end

# struct SingleChoiceExitSolution <: DDPSolution # LKE = LearningKExit
#
#        meshValFun::Array{Float64}
#        meshPolFun::Array{Float64}
#        meshExit::Array{Bool}
#        mV0::Array{Float64}
#        mK0::Array{Float64}
#
# end


# function createsolution(p::TwoChoiceVar, meshValFun::Array{Float64},
#                            meshPolFun1::Array{Float64}, meshPolFun2 ::Array{Float64})
#
#         mV0, mPolFunOne0, mPolFunTwo0 = initialendogstatevars(p, meshValFun)
#
#         TwoChoiceVarolution(meshValFun, meshPolFun1, meshPolFun2,
#                                    mV0, mPolFunOne0, mPolFunTwo0)
# end



# function createsolution(p::LearningKExit, meshValFun::Array{Float64},
#                         meshPolFun::Array{Float64}, meshExit::Array{Float64})
#
#     # nz x nmu matrix of optimal initial capital stock
# 	mK0 = getmK0(p, meshValFun)
#
# 	# interpolator for indirect value V0 (z, mu, P=P0)
# 	itp_V0 = interpolate(p.tStateVectors[1:3], meshValFun[:,:,:,end], Gridded(Linear()))
#
# 	## CONSTRUCT V0 INTERPOLATOR
# 	mV0 = zeros(p.nNodes[2], p.nNodes[3]) # z x mu
# 	for iz = 1 : p.nNodes[2], imu = 1: p.nNodes[3]
# 		# println("iz = $iz, imu = $imu")
# 		# println(mK0[iz, imu])
# 		mV0[iz, imu] = itp_V0[mK0[iz, imu], p.tStateVectors[2][iz],
# 							p.tStateVectors[3][imu]] - (1+p.C0)*mK0[iz, imu]
#
# 		# NOW VALUE DOES NOT INCLUDE INITIAL COST C0
# 	end
# 	# V0 = (K = K0, z, mu, P = P0) indirect utility
#
#     SingleChoiceExitSolution(meshValFun, meshPolFun, meshExit, mK0, mV0)
# end
