using Test
using GBGA
using Evolutionary
using MolecularGraph

# lots of these are by reference to the tests in the Evolutionary API tests
# https://github.com/wildart/Evolutionary.jl/blob/eff5e0a23be5186d1d46bed298c263e8cedd8ff1/test/interface.jl

mol=MolecularGraph.smilestomol("C")

struct TestOptimizer <: Evolutionary.AbstractOptimizer end
method = TestOptimizer()

import Evolutionary.population_size
function population_size(method::TestOptimizer)
    return 1
end

pop=Evolutionary.initial_population(method, mol)
@test length(pop) == Evolutionary.population_size(method)

@test Evolutionary.NonDifferentiable(method, mol)


