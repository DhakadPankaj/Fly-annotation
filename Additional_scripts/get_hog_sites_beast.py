import random
import sys
import subprocess

def subsample_codons(input_file, output_file, num_codons):
    with open(input_file, 'r') as f_in:
        lines = f_in.readlines()

    num_sequences = len(lines) // 2
    num_positions = len(lines[1].strip())

    if num_codons > num_positions // 3:
        print("Error: Number of codons to subsample exceeds the number of available codons. Check if the input fasta file is a single-line fasta file. Exiting...")
        return

    selected_positions = random.sample(range(num_positions // 3), num_codons)
    selected_positions.sort()

    with open(output_file, 'w') as f_out:
        for i in range(num_sequences):
            sequence_name = lines[i * 2].strip()
            sequence = lines[i * 2 + 1].strip()

            selected_codons = [sequence[j * 3 : (j + 1) * 3] for j in selected_positions]
            selected_sequence = ''.join(selected_codons)

            f_out.write(sequence_name + '\n')
            f_out.write(selected_sequence + '\n')

    print("Subsampling of codons completed successfully.")

# Usage example
alignment_file = sys.argv[1]
output_file = sys.argv[2]
num_codons = int(sys.argv[3])

# get sites with less than 5% gaps
#subprocess.run(["clipkit ", alignment_file, "-o ", "supermatrix_filtered.fasta", "-m ", "gappy", "-c ", "-g ", "0.05 ", "-if ", "fasta"], shell=True)

# convert multiline fasta to single line
#subprocess.run(["fasta_formatter ", "-i ", "supermatrix_filtered.fasta", " > ", "single_line.fasta"], shell=True)

# Subsample codons
#input_file = "single_line.fasta"
subsample_codons(alignment_file, output_file, num_codons)