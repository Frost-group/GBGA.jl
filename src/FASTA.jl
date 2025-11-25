# FASTA (protein sequence) genetic algorithm stuff
using Evolutionary
using Random
using StatsBase

export AA_ALPHABET, AA_FREQUENCIES, random_peptide, mutate_peptide, crossover_peptides, peptide_charge, peptide_fitness, peptide_create_ga_config, peptide_optimize

# Amino acid alphabet (20 standard amino acids)
const AA_ALPHABET = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Natural frequency/abundance of amino acids in proteins (in %)
# Based on Swiss-Prot database statistics
const AA_SWISSPROT_FREQUENCIES = Dict(
    'A' => 8.25,   # Alanine
    'C' => 1.37,   # Cysteine
    'D' => 5.45,   # Aspartic acid
    'E' => 6.75,   # Glutamic acid
    'F' => 3.86,   # Phenylalanine
    'G' => 7.07,   # Glycine
    'H' => 2.27,   # Histidine
    'I' => 5.96,   # Isoleucine
    'K' => 5.84,   # Lysine
    'L' => 9.66,   # Leucine
    'M' => 2.42,   # Methionine
    'N' => 4.07,   # Asparagine
    'P' => 4.70,   # Proline
    'Q' => 3.93,   # Glutamine
    'R' => 5.53,   # Arginine
    'S' => 6.56,   # Serine
    'T' => 5.34,   # Threonine
    'V' => 6.87,   # Valine
    'W' => 1.08,   # Tryptophan
    'Y' => 2.92    # Tyrosine
)

# Igor extracted from DBAASP; not sure if short linear cationic only etc.
const AA_DBAASP_FREQUENCIES = Dict(
    'A' => 7.006,   # Alanine
    'C' => 2.535,   # Cysteine
    'D' => 1.266,   # Aspartic acid
    'E' => 1.404,   # Glutamic acid
    'F' => 5.272,   # Phenylalanine
    'G' => 6.976,   # Glycine
    'H' => 1.963,   # Histidine
    'I' => 7.112,   # Isoleucine
    'K' => 16.128,  # Lysine
    'L' => 12.530,  # Leucine
    'M' => 0.929,   # Methionine
    'N' => 2.065,   # Asparagine
    'P' => 4.053,   # Proline
    'Q' => 1.808,   # Glutamine
    'R' => 10.936,  # Arginine
    'S' => 3.410,   # Serine
    'T' => 2.355,   # Threonine
    'V' => 5.540,   # Valine
    'W' => 4.843,   # Tryptophan
    'Y' => 1.868    # Tyrosine
)

function random_peptide(length::Int)
    return join(rand(AA_ALPHABET, length))
end

# Modified random_peptide function to use natural frequencies
function random_peptide_weighted(length::Int; dict=AA_SWISSPROT_FREQUENCIES)
    # Dereference dict into set of weights as vector (ordered as in AA_ALPHABET)
    weights = [dict[aa] for aa in AA_ALPHABET] |> StatsBase.FrequencyWeights
    return sample(AA_ALPHABET, weights, length) |> join
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

# Evolutionary.jl requires zero() for the element type
#  ==> This is some dark Julia magic that I don't fully understand

Base.zero(::Type{Char}) = '\0'
# doing this as a null character, sort of hack around Evolutionary.jl

# Fitness function for a single sequence
function peptide_fitness(sequence::String)
    # qudratic to +5 cation peptide

    fitness = (peptide_charge(sequence) - 5.0)^2
    return fitness
end

# Fitness wrapper for Vector{Char} (required by Evolutionary.jl)
peptide_fitness(seq::Vector{Char}) = peptide_fitness(String(seq))

# Mutation wrapper for Vector{Char} (rng kwarg required by Evolutionary.jl)
function mutate_peptide_vec(seq::Vector{Char}; rng=nothing)
    collect(mutate_peptide(String(seq)))
end

# Crossover wrapper for Vector{Char} (rng kwarg required by Evolutionary.jl)
function crossover_peptides_vec(seq1::Vector{Char}, seq2::Vector{Char}; rng=nothing)
    o1, o2 = crossover_peptides(String(seq1), String(seq2))
    return collect(o1), collect(o2)
end

# Main optimization function
function peptide_optimize(sequence_length::Int;
                         population_size::Int=100,
                         generations::Int=10)
    # Create initial population as Vector{Vector{Char}} for Evolutionary.jl compatibility
    population = [collect(random_peptide(sequence_length)) for _ in 1:population_size]
    
    ga = GA(populationSize = population_size,
            selection = ranklinear(1.5),
            mutation = mutate_peptide_vec,
            crossover = crossover_peptides_vec,
            crossoverRate = 0.8,
            mutationRate = 0.1,
            ε = 0.05)

    result = Evolutionary.optimize(peptide_fitness,
                                 Evolutionary.NoConstraints(),
                                 ga,
                                 population,
                                 Evolutionary.Options(iterations=generations))
   
        @show result

    return String(Evolutionary.minimizer(result))
end

