
using QuantEcon: tauchen

model = :CapitalAdjustModel
dipar_1 = Dict(:nK => 100, :nz => 20, :β => 0.9, :ρ => 0.6, :σ => 0.3)
p_1_all = createmodel(model; dipar_1..., intdim = :All)
p_1_sep = createmodel(model; dipar_1..., intdim = :Separable)
p_1_sep_states = createmodel(model; dipar_1..., intdim = :Separable_States)
p_1_sep_exogstates = createmodel(model; dipar_1..., intdim = :Separable_ExogStates)

T_1_all = transitionmatrix(p_1_all)
T_1_sep = transitionmatrix(p_1_sep)
T_1_sep_states = transitionmatrix(p_1_sep_states)
T_1_sep_exogstates = transitionmatrix(p_1_sep_exogstates)

model = :CapitalAdjustModel2
dipar_2 = Dict(:nK => 20, :nN=>15, :nz => 3, :β => 0.9, :ρ => 0.6, :σ => 0.3)
p_2_all = createmodel(model; dipar_2..., intdim = :All)
p_2_sep = createmodel(model; dipar_2..., intdim = :Separable)
p_2_sep_states = createmodel(model; dipar_2..., intdim = :Separable_States)
p_2_sep_exogstates = createmodel(model; dipar_2..., intdim = :Separable_ExogStates)

T_2_all = transitionmatrix(p_2_all)
T_2_sep = transitionmatrix(p_2_sep)
T_2_sep_states = transitionmatrix(p_2_sep_states)
T_2_sep_exogstates = transitionmatrix(p_2_sep_exogstates)

@testset "transition matrix" begin

	@testset "compare" begin

		@testset "CapitalAdjustModel" begin

			T_all = T_1_all
			T_sep = T_1_sep
			T_sep_states = T_1_sep_states
			T_sep_exogstates = T_1_sep_exogstates

			nK = dipar_1[:nK]
			nz = dipar_1[:nz]

			# compare All and Separable
			# to do

			# compare Separable and Separable_States
			# to do

			# compare Separable_States and Separable_ExogStates
			tT_test = [reshape(T_sep_states,(nK, nz, nz))[j,:,:] for j = 1 : 5]
			@test all([T_sep_exogstates == T_test for T_test in tT_test])

		end

		@testset "CapitalAdjustModel2" begin

			T_all = T_2_all
			T_sep = T_2_sep
			T_sep_states = T_2_sep_states
			T_sep_exogstates = T_2_sep_exogstates

			nK = dipar_2[:nK]
			nN = dipar_2[:nN]
			nz = dipar_2[:nz]

			# compare All and Separable
			# to do

			# compare Separable and Separable_States
			# to do

			# compare Separable_States and Separable_ExogStates
			tT_test = [reshape(T_sep_states, (nK,nN,nz,nz))[j,j,:,:] for j = 1 : 5]
			@test all([T_sep_exogstates == T_test for T_test in tT_test])

		end

	end # compare

	@testset "external transitionmatrix" begin

		nz = dipar_1[:nz]
		ρ = dipar_1[:ρ]
		σ = dipar_1[:σ]

		# have same nodes as in tauchen?
		@test isapprox(collect(tauchen(nz, ρ, σ).state_values),
			p_1_sep_exogstates.tStateVectors[2], rtol = mytol)
		solve(p_1_sep_exogstates, mTransition = tauchen(nz, ρ, σ).p)

	end

end # transition matrix
