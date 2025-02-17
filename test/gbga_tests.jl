using Evolutionary
using MolecularGraph
using Graphs

@testset "GBGA Tests" begin
    # lots of these are by reference to the tests in the Evolutionary API tests
    # https://github.com/wildart/Evolutionary.jl/blob/eff5e0a23be5186d1d46bed298c263e8cedd8ff1/test/interface.jl

    mol = MolecularGraph.smilestomol("C")

    struct TestOptimizer <: Evolutionary.AbstractOptimizer end
    method = TestOptimizer()

    # Define population_size for our TestOptimizer
    Evolutionary.population_size(method::TestOptimizer) = 1

    pop = Evolutionary.initial_population(method, mol)
    @test length(pop) == Evolutionary.population_size(method)

    @test Evolutionary.NonDifferentiable(method, mol)
end 
