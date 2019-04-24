
using QuantEcon: tauchen

model = :NeoClassicalSimple
dipar_neo = Dict(:nK => 100, :nz => 20, :β => 0.9, :ρ => 0.6, :σ => 0.3)
p_neo_all = createmodel(model; dipar_neo..., intdim = :All)
p_neo_sep = createmodel(model; dipar_neo..., intdim = :Separable)
p_neo_sep_states = createmodel(model; dipar_neo..., intdim = :Separable_States)
p_neo_sep_exogstates = createmodel(model; dipar_neo..., intdim = :Separable_ExogStates)

T_neo_all = transitionmatrix(p_neo_all)
T_neo_sep = transitionmatrix(p_neo_sep)
T_neo_sep_states = transitionmatrix(p_neo_sep_states)
T_neo_sep_exogstates = transitionmatrix(p_neo_sep_exogstates)

model = :Intangible
dipar_int = Dict(:nK => 20, :nN=>15, :nz => 3, :β => 0.9, :ρ => 0.6, :σ => 0.3)
p_int_all = createmodel(model; dipar_int..., intdim = :All)
p_int_sep = createmodel(model; dipar_int..., intdim = :Separable)
p_int_sep_states = createmodel(model; dipar_int..., intdim = :Separable_States)
p_int_sep_exogstates = createmodel(model; dipar_int..., intdim = :Separable_ExogStates)

T_int_all = transitionmatrix(p_int_all)
T_int_sep = transitionmatrix(p_int_sep)
T_int_sep_states = transitionmatrix(p_int_sep_states)
T_int_sep_exogstates = transitionmatrix(p_int_sep_exogstates)

@testset "transition matrix" begin

	@testset "compare" begin

		@testset "Neoclassical" begin

			T_all = T_neo_all
			T_sep = T_neo_sep
			T_sep_states = T_neo_sep_states
			T_sep_exogstates = T_neo_sep_exogstates

			nK = dipar_neo[:nK]
			nz = dipar_neo[:nz]

			# compare All and Separable
			# to do

			# compare Separable and Separable_States
			# to do

			# compare Separable_States and Separable_ExogStates
			tT_test = [reshape(T_sep_states,(nK, nz, nz))[j,:,:] for j = 1 : 5]
			@test all([T_sep_exogstates == T_test for T_test in tT_test])

		end

		@testset "Intangible" begin

			T_all = T_int_all
			T_sep = T_int_sep
			T_sep_states = T_int_sep_states
			T_sep_exogstates = T_int_sep_exogstates

			nK = dipar_int[:nK]
			nN = dipar_int[:nN]
			nz = dipar_int[:nz]

			# compare All and Separable
			# to do

			# compare Separable and Separable_States
			# to do

			# compare Separable_States and Separable_ExogStates
			tT_test = [reshape(T_sep_states, (nK,nN,nz,nz))[j,j,:,:] for j = 1 : 5]
			@test all([T_sep_exogstates == T_test for T_test in tT_test])

		end

	end # compare

	@testset "external" begin

		nz = dipar_neo[:nz]
		ρ = dipar_neo[:ρ]
		σ = dipar_neo[:σ]

		# have same nodes as in tauchen?
		@test isapprox(collect(tauchen(nz, ρ, σ).state_values),
			p_neo_sep_exogstates.tStateVectors[2], rtol = mytol)
		solve(p_neo_sep_exogstates, mTransition = tauchen(nz, ρ, σ).p)

	end

end # transition matrix
