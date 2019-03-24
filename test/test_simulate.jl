using DiscreteDynamicProgramming
using Test

# todo: find a better way to test rather than just calling the function

p_neo = createmodel(:NeoClassicalSimple; F = 0.)
optdict = Dict(
	:intdim=>:separable,
	:monotonicity=>true,
	:concavity=>true,
	:rewardmat=>:prebuild_partial,
	:initialize_exact=>true)

shocks_neo = drawshocks(p_neo, 20, 100)
sol_neo = solve(p_neo; optdict...)

@testset "Simulate Neoclassical" begin
	simulate(p_neo, sol_neo, shocks_neo, initialize_exact=false)
	simulate(p_neo, sol_neo, shocks_neo, initialize_exact=true)
end


p_int = createmodel(:Intangible)
optdict = Dict(
	:intdim=>:separable,
	:monotonicity=>true,
	:concavity=>true,
	:rewardmat=>:prebuild_partial,
	:initialize_exact=>true
	)

sol_int = solve(p_int; optdict...)
shocks_int = drawshocks(p_int, 20, 100)

@testset "Simulate Intangible" begin
	simulate(p_int, sol_int, shocks_int, initialize_exact=false)
	simulate(p_int, sol_int, shocks_int, initialize_exact=true)
end
