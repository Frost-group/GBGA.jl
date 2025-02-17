

@testset "FASTA Genetic Algorithm Tests" begin
    @testset "Random Peptide Generation" begin
        # Test length
        seq = random_peptide(10)
        @test length(seq) == 10
        
        # Test valid amino acids
        @test all(aa -> aa ∈ AA_ALPHABET, seq)
        
        # Test randomness (different sequences)
        seq2 = random_peptide(10)
        @test seq != seq2  # Nb: fails once in 22^10 blue moons
    end

    @testset "Mutation Operator" begin
        seq = "ACDEFGHIKL"
        mutated = mutate_peptide(seq, 1.0)  # 100% mutation rate
        
        # Test length preservation
        @test length(mutated) == length(seq)
        
        # Test valid amino acids after mutation
        @test all(aa -> aa ∈ AA_ALPHABET, mutated)
        
        # Test no mutation case
        no_mutation = mutate_peptide(seq, 0.0)
        @test no_mutation == seq
    end

    @testset "Crossover Operator" begin
        seq1 = "ACDEFGHIKL"
        seq2 = "MNPQRSTVWY"
        
        offspring1, offspring2 = crossover_peptides(seq1, seq2)
        
        # Test length preservation
        @test length(offspring1) == length(seq1)
        @test length(offspring2) == length(seq2)
        
        # Test valid amino acids after crossover
        @test all(aa -> aa ∈ AA_ALPHABET, offspring1)
        @test all(aa -> aa ∈ AA_ALPHABET, offspring2)
        
        # Test different length sequences throw error
        @test_throws AssertionError crossover_peptides("ABC", "ABCD")
    end

    @testset "Peptide Charge Calculation" begin
        # Test basic amino acids
        @test peptide_charge("R") ≈ peptide_charge("K") ≈ 1.0  # Lysine
        @test peptide_charge("H") ≈ 0.1  # Histidine
        
        # Test acidic amino acids
        @test peptide_charge("D") ≈ -1.0  # Aspartic acid
        @test peptide_charge("E") ≈ -1.0  # Glutamic acid
        
        # Test neutral amino acids
        @test peptide_charge("A") ≈ 0.0  # Alanine
        
        # Test combinations
        @test peptide_charge("RK") ≈ 2.0    # Two positive
        @test peptide_charge("DE") ≈ -2.0   # Two negative
        @test peptide_charge("RD") ≈ 0.0    # One positive, one negative
        @test peptide_charge("RKDE") ≈ 0.0  # Two positive, two negative

        @test peptide_charge("DEADTEEF") ≈ -5.0
    end



    @testset "Optimization" begin
        # Basic test to ensure optimization runs
        result = peptide_optimize(5, generations=10)
        
        # Test result properties
        @test length(result) == 5
        @test all(aa -> aa ∈ AA_ALPHABET, result)
    end
end 