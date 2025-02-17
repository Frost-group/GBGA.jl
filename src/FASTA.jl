# FASTA (protein sequence) genetic algorithm stuff
using Evolutionary
using Random

export AA_ALPHABET, random_peptide, mutate_peptide, crossover_peptides, peptide_charge, peptide_fitness, peptide_create_ga_config, peptide_optimize

# Amino acid alphabet (20 standard amino acids)
const AA_ALPHABET = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Generate a random peptide sequence of given length
function random_peptide(length::Int)
    return join(rand(AA_ALPHABET, length))
end

# Mutation operator: randomly select new amino acid for each position with probability mutation_rate
function mutate_peptide(sequence::String, mutation_rate::Float64=0.01)
    chars = collect(sequence)
    for i in 1:length(chars)
        if rand() < mutation_rate
            chars[i] = rand(AA_ALPHABET)
        end
    end
    return join(chars)
end

# Crossover operator: single-point crossover between two peptide sequences
function crossover_peptides(seq1::String, seq2::String)
    @assert length(seq1) == length(seq2) "Sequences must be of equal length"
    point = rand(1:length(seq1))
    offspring1 = seq1[1:point] * seq2[point+1:end]
    offspring2 = seq2[1:point] * seq1[point+1:end]
    return offspring1, offspring2
end

# total charge
function peptide_charge(sequence::String)
    # Amino acid charges at physiological pH (7.4)
    # Positive: R (Arg), K (Lys), H (His)
    # Negative: D (Asp), E (Glu)
    charges = Dict(
        'R' => 1.0,  # Arginine
        'K' => 1.0,  # Lysine
        'H' => 0.1,  # Histidine (partially charged at pH 7.4)
        'D' => -1.0, # Aspartic acid
        'E' => -1.0, # Glutamic acid
        # All other amino acids are neutral (charge = 0)
    )
    
    # Calculate total charge
    total_charge = sum(get(charges, aa, 0.0) for aa in sequence)
    return total_charge
end

#  Evolutionary.jl requires a zero function for the population type
#  ==> This is some dark Julia magic that I don't fully understand
# Add Base.zero for String type to satisfy Evolutionary.jl requirements
Base.zero(::Type{String}) = ""
# Add Base.zero for Vector{Float64} type
Base.zero(::Type{Vector{Float64}}) = Float64[]

# Vector version of fitness function; sort of a DIY dot operator
function peptide_fitness(population::Vector{String})
    println("population: $population")
    return Float64[peptide_fitness(seq) for seq in population]
end

# per sequence version of fitness function; where the complexity goes
function peptide_fitness(sequence::String)
    return peptide_charge(sequence)
end

# Main optimization function
function peptide_optimize(sequence_length::Int;
                         population_size::Int=100,
                         generations::Int=100)
    # Create initial population
    population = [random_peptide(sequence_length) for _ in 1:population_size]
    
    # Define GA parameters
    ga = GA(populationSize = population_size,
          #  selection = tournament(3),
            mutation = x -> mutate_peptide(x),
            crossover = (x, y) -> crossover_peptides(x, y)[1],
            crossoverRate = 0.8,
            mutationRate = 0.1,
            ε = 0.05)

    result = Evolutionary.optimize(peptide_fitness,
                                 population,
                                 ga,
                                 Evolutionary.Options(iterations=generations))
    
    return Evolutionary.minimizer(result)
end

