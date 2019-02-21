
using DiscreteDynamicProgramming
include("helpfunctions.jl")
using Test
mytol = 1e-4

# CHANGE TO F > 0

model = :NeoClassicalSimple
dipar = Dict(:nK => 40, :nz => 11, :γ => 0.5, :F => 0., :τ => 0.3,
    :intdim => :separable)
p_neoclassical = createmodel(model; dipar...,
    monotonicity = false, concavity = false)
p_neoclassical_mon = createmodel(model; dipar...,
    monotonicity = true, concavity = false)
p_neoclassical_conc = createmodel(model; dipar...,
    monotonicity = false, concavity = true)

sol_neoclassical = solve(p_neoclassical)
sol_neoclassical_mon = solve(p_neoclassical_mon)
sol_neoclassical_conc = solve(p_neoclassical_conc)

@testset "monotonicity concavity" begin

    compare_solutions("Neoclassical_monotonicity", sol_neoclassical, sol_neoclassical_mon)
    compare_solutions("Neoclassical_concavity", sol_neoclassical, sol_neoclassical_conc)

end
