module GBGA

using MolecularGraph
using Evolutionary
using Graphs
import OpenSMILES 
import Evolutionary: population_size  # Import population_size

greet() = print("Hello World!")

# Override Evolutionary initial population
import Evolutionary.initial_population
function initial_population(method::M, individual::MolGraph) where {M<:Evolutionary.AbstractOptimizer}
    # Return an array of molecules (population) of the specified size
    return [MolecularGraph.smilestomol("Clc4cc2c(C(/c1ncccc1CC2)=C3/CCNCC3)cc4") 
            for _ in 1:population_size(method)]
end

import Evolutionary.NonDifferentiable
NonDifferentiable(f, x::MolGraph) = NonDifferentiable{Real,typeof(x)}(f, f(x), deepcopy(x), [0,])

end # module
