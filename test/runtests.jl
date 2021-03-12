using Test
using GBGA
using Evolutionary
using MolecularGraph

# lots of these are by reference to the tests in the Evolutionary API tests
# https://github.com/wildart/Evolutionary.jl/blob/eff5e0a23be5186d1d46bed298c263e8cedd8ff1/test/interface.jl

m=MolecularGraph.smilestomol("C")

struct TestOptimizer <: Evolutionary.AbstractOptimizer end
mthd = TestOptimizer()

pop=Evolutionary.initial_population(mthd, m)
@test length(pop) == Evolutionary.population_size(mthd)

@test Evolutionary.NonDifferentiable(mthd, m)


