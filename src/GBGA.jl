module GBGA

import OpenSMILES 
import Evolutionary

greet() = print("Hello World!")

# Override Evolutionary initial population
import Evolutionary.initial_population
function initial_population(method::M, individual::SimpleWeightedGraph) where {M<:Evolutionary.AbstractOptimizer}
    OpenSMILES.ParseSMILES("Clc4cc2c(C(/c1ncccc1CC2)=C3/CCNCC3)cc4")
end

import Evolutionary.NonDifferentiable
NonDifferentiable(f, x::GraphMol) = NonDifferentiable{Real,typeof(x)}(f, f(x), deepcopy(x), [0,])


end # module
