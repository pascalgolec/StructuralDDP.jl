
using DiscreteDynamicProgramming
using Test
mytol = 1e-4

model = :NeoClassicalSimple
nK = 100
nz = 20
dipar = Dict(:nK => nK, :nz => nz, :β => 0.9, :ρ => 0.6, :σ => 0.3)
p_neo_all = createmodel(model; dipar..., intdim = :All)
p_neo_sep = createmodel(model; dipar..., intdim = :Separable)
p_neo_sep_states = createmodel(model; dipar..., intdim = :Separable_States)
p_neo_sep_exogstates = createmodel(model; dipar..., intdim = :Separable_ExogStates)

T_neo_all = transitionmatrix(p_neo_all)
T_neo_sep = transitionmatrix(p_neo_sep)
T_neo_sep_states = transitionmatrix(p_neo_sep_states)
T_neo_sep_exogstates = transitionmatrix(p_neo_sep_exogstates)

#
# p_intangible_sep = createmodel(NeoClassicalIntangible, intdim = :separable,
# 		β  = 0.9, ρ  = 0.6, σ  = 0.3, nK = 20, nN = 15, nz = 3, nShocks = 3)
#
# p_intangible_int = createmodel(NeoClassicalIntangible, intdim = :intermediate,
# 		β  = 0.9, ρ  = 0.6, σ  = 0.3, nK = 20, nN = 15, nz = 3, nShocks = 3)
#
# p_intangible_SA = createmodel(NeoClassicalIntangible, intdim = :SA,
# 		β  = 0.9, ρ  = 0.6, σ  = 0.3, nK = 20, nN = 15, nz = 3, nShocks = 3)
#
# T_intangible_sep = transitionmatrix(p_intangible_sep)
# T_intangible_int = transitionmatrix(p_intangible_int)
# T_intangible_SA = transitionmatrix(p_intangible_SA)
#
#
#
# p_lev_sep = createmodel(FinancingLeverageSimple, intdim = :separable,
#         β  = 0.9, ρ = 0.6, σ = 0.3, nK = 50, nLev = 7, nz = 5, nShocks = 5)
#
# p_lev_int = createmodel(FinancingLeverageSimple, intdim = :intermediate,
#         β  = 0.9, ρ = 0.6, σ = 0.3, nK = 50, nLev = 7, nz = 5, nShocks = 5)
#
# p_lev_SA = createmodel(FinancingLeverageSimple, intdim = :SA,
#         β  = 0.9, ρ = 0.6, σ = 0.3, nK = 50, nLev = 7, nz = 5, nShocks = 5)
#
# T_lev_sep = transitionmatrix(p_lev_sep)
# T_lev_int = transitionmatrix(p_lev_int)
# T_lev_SA = transitionmatrix(p_lev_SA)
#
#
# p_IEI_int = createmodel(IEI, intdim = :intermediate,
#         β  = 0.9, ρ = 0.6, σ = 0.3, nK = 10, nS = 7, nz = 5, nShocks = [3,10])
#
# p_IEI_SA = createmodel(IEI, intdim = :SA,
#         β  = 0.9, ρ = 0.6, σ = 0.3, nK = 10, nS = 7, nz = 5, nShocks = [3,10])
#
# T_IEI_int = transitionmatrix(p_IEI_int)
# T_IEI_SA = transitionmatrix(p_IEI_SA)


@testset "transition matrix" begin

	# @testset "separable" begin
	#
	# 	@testset "Neoclassical" begin
	# 		@test T_neoclassical_sep[1,1] ≈ 0.16666666666666674 rtol=mytol #
	# 		@test T_neoclassical_sep[2,1] ≈ 0.16464367430685903 rtol=mytol #
	# 		@test T_neoclassical_sep[4,1] ≈ 0.0 rtol=mytol #
	# 		@test T_neoclassical_sep[4,3] ≈ 0.03535632569314069 rtol=mytol
	# 		@test T_neoclassical_sep[1,5] ≈ 0.5333333333333334 rtol=mytol #
	# 	end
	#
	# 	@testset "Intangible" begin
	# 		@test T_intangible_sep[1,1] ≈ 0.5896866307747166 rtol=mytol #
	# 		@test T_intangible_sep[2,1] ≈ 0.07698003589195013 rtol=mytol #
	# 		@test T_intangible_sep[3,2] ≈ 0.4103133692252835 rtol=mytol #
	# 	end
	#
    #     @testset "Financing" begin
	# 		@test T_lev_sep[1,1] ≈ 0.32290029376229984 rtol=mytol #
	# 		@test T_lev_sep[2,1] ≈ 0.08298807389226581 rtol=mytol #
	# 		@test T_lev_sep[4,1] ≈ 0.0 rtol=mytol #
	# 		@test T_lev_sep[4,3] ≈ 0.40895252244807356 rtol=mytol
	# 		@test T_lev_sep[1,5] ≈ 0.0 rtol=mytol #
	# 	end
	# end


	@testset "compare" begin

		@testset "Neoclassical" begin

			T_all = T_neo_all
			T_sep = T_neo_sep
			T_sep_states = T_neo_sep_states
			T_sep_exogstates = T_neo_sep_exogstates

			# compare All and Separable
			# to do

			# compare Separable and Separable_States
			# to do

			# compare Separable_States and Separable_ExogStates
			tT_test = [reshape(T_sep_states,(nK, nz, nz))[j,:,:] for j = 1 : 5]
			@test all([T_sep_exogstates == T_test for T_test in tT_test])

		end

		# @testset "Intangible" begin
        # 	T_test = reshape(T_intangible_int,(20, 15, 3, 3))[2, 4, :,:]
		# 	@test all(T_intangible_sep == T_test)
		# 	T_test = reshape(T_intangible_int,(20, 15, 3, 3))[10, 5, :,:]
		# 	@test all(T_intangible_sep == T_test)
		#
		# 	# contract SA to get intermediate
		# 	test = T_intangible_SA[1:20*15*3,:]
		# 	test = Matrix(reshape(test, (20*15*3, 20, 15, 3))[:,1,1,:])
		# 	@test test  ≈ Matrix(T_intangible_int) rtol=mytol
		# end
		#
        # @testset "Financing" begin
		# 	# contract intermediate to get separate
        #     T_test = reshape(T_lev_int,(50, 7, 5, 5))[2, 4, :,:]
		# 	@test all(T_lev_sep == T_test)
		# 	T_test = reshape(T_lev_int,(50, 7, 5, 5))[10, 5, :,:]
		# 	@test all(T_lev_sep == T_test)
		#
		# 	# contract SA to get intermediate
		# 	test = T_lev_SA[1:50*7*5,:]
		# 	test = Matrix(reshape(test, (50*7*5, 50, 7, 5))[:,1,1,:])
		# 	@test test  ≈ Matrix(T_lev_int) rtol=mytol
		# end
		#
		# @testset "IEI" begin
		# 	# contract SA to get intermediate
		# 	test = T_IEI_SA[1:10*7*5,:]
		# 	test = Matrix(reshape(test, (10*7*5, 10, 7, 5))[:,1,1,:])
		# 	@test test  ≈ Matrix(T_IEI_int) rtol=mytol
		# end

	end

end # transition matrix
