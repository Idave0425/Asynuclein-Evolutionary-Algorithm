import random
import subprocess
import tempfile
import sys
import os
import csv

POSSIBLE_AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
ALPHA_SY = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
# Population Initialization
# NAC region has 35 AA, we choose ~10 before and after to make 55 AA molecule
# Call this function like initialize_population(30, 55)
# def initialize_population(n, m):
#  pop_i = [[random.choi(POSSIBLE_AA) for _ in range(m)] for _ in range(cen)]
#  return pop_i

def initialize_population(n, m, alpha_seq):
    """
    Generate a population of mutated Alpha-Synuclein sequences.

    Args:
    - n (float): Percentage of amino acids to mutate (0-100).
    - m (int): Length of each sequence (e.g., 55 AA).
    - alpha_seq (str): Amino acid sequence of Alpha-Synuclein.

    Returns:
    - list: A list of mutated sequences.
    """
    population_size = m  # Number of mutated sequences to generate
    num_mutations = int(len(alpha_seq) * (n / 100))  # Calculate the number of mutations
    population = []

    for _ in range(population_size):
        mutated_seq = list(alpha_seq)  # Convert sequence to a list for mutability
        mutation_indices = random.sample(range(len(alpha_seq)), num_mutations)  # Select mutation positions
        for index in mutation_indices:
            # Ensure the new amino acid is different from the original
            new_aa = random.choice([aa for aa in POSSIBLE_AA if aa != mutated_seq[index]])
            mutated_seq[index] = new_aa
        population.append("".join(mutated_seq))  # Add mutated sequence to the population

    return population

def get_likelihood_of_interaction(tsv_file_path):
    with open(tsv_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        
        # Read through the file rows and capture LikelihoodOfInteraction value from the last row
        likelihood_of_interaction = None
        for row in reader:
            likelihood_of_interaction = float(row['LikelihoodOfInteraction'])
        
        return likelihood_of_interaction

def generate_mutated_fasta(original_fasta_path, mutated_sequence, output_fasta_path='mutated.fasta'):
    try:
        with open(original_fasta_path, 'r') as infile, open(output_fasta_path, 'w') as outfile:
            species = None
            sequence = ""

            for line in infile:
                line = line.strip()

                # Check for a new species identifier
                if line.startswith('>'):
                    # If there's a previous species, write its sequence
                    if species:
                        outfile.write(f"{species}\n")
                        for i in range(0, len(sequence), 60):
                            outfile.write(f"{sequence[i:i+60]}\n")

                    # Start a new species
                    species = line
                    sequence = ""

                else:
                    sequence += line

            # Handle the last species
            if species:
                if 'Homo_sapiens' in species:
                    outfile.write(f"{species}\n")
                    for i in range(0, len(mutated_sequence), 60):
                        outfile.write(f"{mutated_sequence[i:i+60]}\n")
                else:
                    outfile.write(f"{species}\n")
                    for i in range(0, len(sequence), 60):
                        outfile.write(f"{sequence[i:i+60]}\n")

        print(f"Fasta file successfully saved as {output_fasta_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

def calculate_specificity(individual):
    # Generate Mutated Fasta
    aa_sequence_str = ''.join(individual)
    generate_mutated_fasta('/Users/ishandave/Documents/asynuclein_evoalgo/original_asyn.fasta', aa_sequence_str)
    # Build Pic Command
    pic_command = "python ./pic_specificity_algorithm/mi_smart_filters.py -f1 /Users/ishandave/Documents/asynuclein_evoalgo/original_asyn.fasta -f2 /Users/ishandave/Documents/asynuclein_evoalgo/mutated.fasta -p1 SNCA -p2 SNCA_M -g v -m Homo_sapiens -s /Users/ishandave/Documents/asynuclein_evoalgo/interactingProteins.tsv"
    ## Call Pic command
    os.system(pic_command)
    specificity_score = get_likelihood_of_interaction('interactingProteins.tsv')
    return specificity_score

# Function to remove x individuals with lowest scores
def highestscore_selection(population, scores, x):
    sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    updated_population = [population[i] for i in sorted_indices[:x]]
    return updated_population


# Function to generate new individuals near existing ones
def generate_new_individuals(population, pop_size, X):
    new_generation = population[:]
    while len(new_generation) < pop_size:
        parent = random.choice(population)
        new_individual = list(parent)
        indices_to_change = random.sample(range(len(parent)), X)
        for index in indices_to_change:
            new_individual[index] = random.choice(POSSIBLE_AA)  # Assuming POSSIBLE_AA is a list of amino acids
        new_generation.append(new_individual)  # Append the list of amino acids directly
    return new_generation

# Iterative Filtering and Modification
def optimize_molecule(num_gen, pop_size, aa_len, sel_str, mut_rate):
    pop_i = initialize_population(50, pop_size, ALPHA_SY)
    pop = pop_i[:]
    likelihoods = []
    for i in range(num_gen):
        print(f"Generation {i}")
        print(f"population {pop}")
        pop_specificities = []
        for individual in pop:
            print(individual)
            specificity_score = calculate_specificity(individual)
            print(specificity_score)
            pop_specificities.append(specificity_score)
        selected_indivuals = highestscore_selection(pop_i, pop_specificities, sel_str)
        new_generation = generate_new_individuals(selected_indivuals, pop_size, mut_rate)
        pop = new_generation[:]
        likelihoods.append(pop_specificities)
        print(pop_specificities)
    return pop, likelihoods


optimize_molecule(1, 5, 55, 1, 130)

