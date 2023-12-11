"""
# Create a Python virtual environment named 'venv'
python3.10 -m venv venv

# Activate the 'venv' virtual environment
source venv/bin/activate

# Install Python packages
pip install -r requirements.txt
"""

import random
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt

BASES: List[str] = ['A', 'C', 'G', 'T']  # List of DNA bases.
ERROR_RATES: List[float] = [0.02, 0.05, 0.10]  # List of error rates to simulate different types of errors in DNA sequences.
ERROR_TYPES: List[str] = ['insertion', 'deletion', 'mismatch']  # Define the types of errors to analyze.
GC_CONTENT: float = 0.6  # GC content of DNA sequences. GC content is the percentage of Guanine (G) and Cytosine (C) bases.
HOMOPOL_BIAS: float = 0.7  # Homopolymer bias represents the bias towards consecutive repetitions of the same base in a sequence.
NUM_SEQUENCES: int = 100  # Number of DNA sequences to generate.
RANDOM_SEED: int = 42  # Random seed for reproducibility of random number generation.
SEQ_LENGTH: int = 100  # Length of each DNA sequence.


def create_random_dna_sequence(length: int = 100, gc_content: float = 0.60) -> str:
    """Generate a random DNA sequence of a given length with a specified GC content.

    :param length: Length of the DNA sequence.
    :param gc_content: The proportion of G and C in the DNA sequence.
    :return: A string representing the DNA sequence.
    """
    gc_count = int(length * gc_content)  # Calculate the number of GC base pairs.
    at_count = length - gc_count  # Calculate the number of AT base pairs.
    bases = ['G', 'C'] * gc_count + ['A', 'T'] * at_count  # Create sequence with balanced GC content
    random.shuffle(bases)  # Shuffle the sequence randomly.

    return ''.join(bases[:length])  # Convert the list of characters into a string and return it.


def introduce_errors_with_tracking(sequence: str, error_rate: float) -> Tuple[str, List[Optional[Tuple[int, str]]]]:
    """Introduces errors into a DNA sequence with advanced bias and tracks the errors.

    :param sequence: The original DNA template sequence.
    :param error_rate: The rate at which errors are introduced.
    :return: Modified sequence with errors and a list of error data (position, type).
    """
    # Identify homopolymer locations
    analysis = analyze_homopolymers(sequences=[sequence])
    homo_positions = [positions for values in analysis.values() for pos_sets in values for positions in pos_sets]

    sequence_length = len(sequence)  # Get the length of the sequence.
    num_errors = int(sequence_length * error_rate)  # Calculate the number of errors based on the error rate.
    error_positions = []
    for error_n in range(num_errors):
        # Generate error positions with weight towards the end of the sequence.
        position = round(random.triangular(low=1, high=sequence_length, mode=sequence_length))
        error_positions.append(position)

        # Check if the position is part of a homopolymer of at least length 4.
        if position in homo_positions and random.random() <= HOMOPOL_BIAS:
            # Generate additional random homopolymer error
            i = homo_positions.index(position)
            if i < len(homo_positions)-1:
                position = homo_positions[i+1]
            else:
                position = homo_positions[i-1]
            error_positions.append(position)
    error_positions.sort()
    error_types = random.choices(population=ERROR_TYPES, k=len(error_positions))

    sequence_nucs = list(sequence)  # Convert the input sequence into a list of characters.
    shift = 0  # Increase or decrease the sequence length
    error_data = []  # Initialize an empty list to track error data.
    for position, error_type in zip(error_positions, error_types):
        i = position - 1 + shift
        if error_type == 'insertion':  # If the error type is insertion:
            sequence_nucs.insert(i, random.choice(BASES))  # Insert a random base at the chosen position.
            shift += 1  # Increase the sequence length.
        elif error_type == 'deletion':  # If the error type is deletion:
            del sequence_nucs[i]  # Delete the base at the chosen position.
            shift -= 1  # Decrease the sequence length.
        elif error_type == 'mismatch':  # If the error type is mismatch:
            # Replace the base at the chosen position with a random non-matching base.
            sequence_nucs[i] = random.choice([b for b in BASES if b != sequence_nucs[i]])
        error_data.append((position, error_type))  # Track the mismatch error.

    # Return converted modified sequence back to a string and return it along with the error data.
    return ''.join(sequence_nucs), error_data


def generate_sequences_with_error_tracking(template_sequence: str,
                                           error_rate: float,
                                           num_sequences=100) -> Dict[str, List[Tuple[int, str]]]:
    """Generates sets of sequences with specified error rates and tracks the errors.

    :param template_sequence: The template DNA sequence.
    :param error_rate: Error rates.
    :param num_sequences: Number of sequences to generate.
    :return: A dictionary with modified sequences as keys and lists of error data as values.
    """
    all_error_data = {}  # Initialize an empty dictionary to store sequences and their error data.
    for _ in range(num_sequences):
        # Generate a modified sequence and track errors.
        modified_sequence, sequence_errors = introduce_errors_with_tracking(sequence=template_sequence,
                                                                            error_rate=error_rate)

        # Store the error data for the current error rate in the dictionary.
        all_error_data[modified_sequence] = sequence_errors

    return all_error_data


def plot_error_distribution(all_error_data):
    error_rates = all_error_data.keys()  # Get the list of error rates from the provided error data.
    colors = {0.02: 'blue', 0.05: 'red', 0.10: 'green'}  # Assign colors to error rates.

    fig, axes = plt.subplots(len(ERROR_TYPES), figsize=(12, 8), sharex=True)  # Create subplots for each error type.

    for i, error_type in enumerate(ERROR_TYPES):  # Iterate over error types.
        ax = axes[i]  # Select the current subplot.
        ax.set_title(f'{error_type.capitalize()} Error Position Distribution')  # Set the title of the subplot.
        ax.set_xlabel('Position in Sequence')  # Set the x-axis label.
        ax.set_ylabel('Frequency of Errors')  # Set the y-axis label.

        # Prepare data for stacked histograms
        histograms = []
        labels = []
        for rate in error_rates:  # Iterate over error rates.
            # Get error data for the current error rate.
            all_errors = [errors for values in all_error_data[rate].values() for errors in values]

            # Extract positions of errors of the current type.
            error_positions = [pos for pos, typ in all_errors if typ == error_type]
            histograms.append(error_positions)
            labels.append(f'{error_type.capitalize()} at {rate * 100}%')

        # Plot stacked histograms
        ax.hist(histograms, bins=range(SEQ_LENGTH+1), alpha=0.7, label=labels, color=[colors[rate] for rate in error_rates], stacked=True)

        ax.legend()  # Add a legend to the subplot.

    plt.tight_layout()  # Adjust subplot layout.
    plt.show()  # Display the plot.


def calculate_gc_content(sequence: str) -> float:
    """Calculates the GC content of a DNA sequence.

    :param sequence: The DNA sequence.
    :return: GC content as a percentage of the total sequence length.
    """
    # Calculate and return the GC content as a percentage.
    return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100


def plot_length_per_error_rate(all_error_data: dict) -> None:
    """Plots the distribution of sequence lengths for each error rate, stacked.

    :param all_error_data: Dictionary with error rates as keys and dictionary of lists of error data as values.
    """
    plt.figure(figsize=(12, 6))  # Create a figure for plotting.

    graph_data = {}
    all_lengths = []
    error_rates = list(all_error_data.keys())  # Get the list of error rates from the provided error data.
    for rate in error_rates:  # Iterate over error rates and corresponding sequence lengths.
        error_data = all_error_data[rate]
        lengths = [len(seq) for seq in error_data.keys()]
        graph_data[rate] = lengths
        all_lengths.extend(lengths)

    min_length = min(all_lengths)  # Find the minimum sequence length.
    max_length = max(all_lengths)  # Find the maximum sequence length.
    bins = range(min_length, max_length + 1)  # Define bins for the histogram.

    # Prepare data for stacked histograms
    histograms = [graph_data[rate] for rate in error_rates]
    labels = [f'{rate * 100}% Error Rate' for rate in error_rates]

    # Plot stacked histograms
    plt.hist(histograms, bins=bins, stacked=True, label=labels)

    plt.title('Sequence Length Distribution per Error Rate (Stacked)')  # Set the title for the plot.
    plt.xlabel('Sequence Length')  # Set the x-axis label.
    plt.ylabel('Frequency')  # Set the y-axis label.
    plt.legend()  # Add a legend.
    plt.show()  # Display the plot.


def calculate_gc_content(sequence: str) -> float:
    """Calculate the GC content of a DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return 100 * gc_count / len(sequence)


def plot_gc_per_error_rate(all_error_data: dict) -> None:
    """Plots the distribution of GC contents for each error rate, stacked.

    :param all_error_data: Dictionary with error rates as keys and dictionary of lists of error data as values.
    """
    plt.figure(figsize=(12, 6))  # Create a figure for plotting.

    graph_data = {}
    all_gc_content = []
    error_rates = list(all_error_data.keys())  # Get the list of error rates from the provided error data.
    for rate in error_rates:  # Iterate over error rates and corresponding gc contents.
        error_data = all_error_data[rate]
        gc_contents = [calculate_gc_content(seq) for seq in error_data.keys()]
        graph_data[rate] = gc_contents
        all_gc_content.extend(gc_contents)

    min_gc_content = min(all_gc_content)  # Find the minimum gc content.
    max_gc_content = max(all_gc_content)  # Find the maximum gc content.
    bins = range(int(min_gc_content), int(max_gc_content) + 2)  # Define bins for the histogram.

    # Prepare data for stacked histograms
    histograms = [graph_data[rate] for rate in error_rates]
    labels = [f'{rate * 100}% Error Rate' for rate in error_rates]

    # Plot stacked histograms
    plt.hist(histograms, bins=bins, stacked=True, label=labels)

    plt.title('GC Content Distribution per Error Rate (Stacked)')  # Set the title for the plot.
    plt.xlabel('GC Content (%)')  # Set the x-axis label.
    plt.ylabel('Frequency')  # Set the y-axis label.
    plt.legend()  # Add a legend.
    plt.show()  # Display the plot.


def analyze_homopolymers(sequences: List[str]) -> Dict[str, List[Optional[List[int]]]]:
    """Analyzes homopolymers in a set of sequences.

    :param sequences: A list of DNA sequences.
    :return dict: A dictionary containing homopolymer analysis data for each nucleotide type.
    """
    # Initialize a dictionary to store homopolymer analysis data for each nucleotide type.
    analysis = {nuc: [] for nuc in BASES}

    # Iterate over the provided sequences.
    for seq in sequences:
        current_homopolymer = None  # Initialize a variable to track the current homopolymer.

        # Iterate over the positions in the sequence.
        for i in range(len(seq)):
            if i == 0 or seq[i] != seq[i - 1]:  # Check if the current position starts a new homopolymer.
                current_homopolymer = (i, seq[i])  # Start tracking a new homopolymer.

            if i == len(seq) - 1 or seq[i] != seq[i + 1]:  # Check if the current position ends the current homopolymer.
                if i - current_homopolymer[0] >= 3:  # Check if the homopolymer length is at least 4.
                    # Append the homopolymer positions to the corresponding nucleotide in the analysis dictionary.
                    analysis[current_homopolymer[1]].append(list(range(current_homopolymer[0], i + 1)))

    # Return the homopolymer analysis data.
    return analysis


def plot_homopolymer_ratios(all_error_data: dict) -> None:
    """
    Plots the homopolymer ratios relative to the sequence length per nucleotide per error rate.

    :param all_error_data: Dictionary with error rates as keys and lists of sequences as values.
    """
    plt.figure(figsize=(12, 6))  # Create a figure for plotting.

    # Initialize a dictionary to hold total homopolymer length and total sequence length for each nucleotide at each error rate
    homopolymer_lengths = {rate: {nuc: 0 for nuc in BASES} for rate in ERROR_RATES}
    total_sequence_lengths = {rate: 0 for rate in ERROR_RATES}

    # Iterate over error rates and sequences to calculate total homopolymer lengths
    for rate, sequences in all_error_data.items():
        for seq in sequences:
            homopolymer_data = analyze_homopolymers([seq])
            seq_length = len(seq)
            total_sequence_lengths[rate] += seq_length

            for nucleotide in BASES:
                # Calculate the total length of all homopolymers for the nucleotide
                homopolymer_lengths[rate][nucleotide] += sum(len(homopoly) for homopoly in homopolymer_data[nucleotide])

    # Plot the homopolymer ratio for each nucleotide at each error rate
    for nucleotide in BASES:
        ratios = [homopolymer_lengths[rate][nucleotide] / total_sequence_lengths[rate] if total_sequence_lengths[rate] > 0 else 0 for rate in ERROR_RATES]
        plt.plot(ERROR_RATES, ratios, label=f'{nucleotide} Homopolymer', marker='o')

    plt.title('Homopolymer Ratios Relative to Sequence Length per Nucleotide per Error Rate')  # Set the title of the plot.
    plt.xlabel('Error Rate')  # Set the x-axis label.
    plt.ylabel('Homopolymer Ratio')  # Set the y-axis label.
    plt.legend()  # Add a legend.
    plt.show()  # Display the plot.


def plot_homopolymer_positions(all_error_data: dict) -> None:
    """Plots the homopolymer positions per error rate per nucleotide in a single figure.

    :param all_error_data: Dictionary with error rates as keys and dictionary of lists of error data as values.
    """
    plt.figure(figsize=(12, 6))  # Create a figure for plotting.

    graph_data = {}
    error_rates = all_error_data.keys()  # Get the list of error rates from the provided error data.
    for rate in error_rates:  # Iterate over error rates.
        error_data = all_error_data[rate]
        sequences = [seq for seq in error_data.keys()]
        graph_data[rate] = analyze_homopolymers(sequences=sequences)

    # Define colors for each nucleotide type.
    nucleotide_colors = {'A': 'red', 'C': 'green', 'G': 'blue', 'T': 'purple'}

    # Create an empty dictionary to store legend handles for each nucleotide.
    legend_handles = {nucleotide: None for nucleotide in nucleotide_colors.keys()}

    for rate in error_rates:  # Iterate over error rates.
        for nucleotide, color in nucleotide_colors.items():  # Iterate over nucleotide types and their colors.
            # Extract homopolymer positions for the current error rate and nucleotide.
            positions = [pos for pos_list in graph_data[rate][nucleotide] for pos in pos_list]

            # Scatter plot of homopolymer positions.
            scatter = plt.scatter(positions, [rate] * len(positions), alpha=0.5, color=color, label=nucleotide)
            
            # Update the legend handle for each nucleotide.
            if legend_handles[nucleotide] is None:
                legend_handles[nucleotide] = scatter

    plt.title('Homopolymer Positions per Error Rate per Nucleotide')  # Set the title of the plot.
    plt.xlabel('Position in Sequence')  # Set the x-axis label.
    plt.ylabel('Error Rate')  # Set the y-axis label.

    # Create a legend using the updated handles and place it between error rates 0.06 to 0.09.
    plt.legend(handles=list(legend_handles.values()), loc='center left', bbox_to_anchor=(1, 0.5), title='Nucleotides')

    plt.grid(True)  # Add grid lines.
    plt.show()  # Display the plot.


def main():
    random.seed(RANDOM_SEED)  # Seed for reproducibility

    # Task1: Create a random DNA sequence
    dna_sequence = create_random_dna_sequence(length=SEQ_LENGTH,
                                              gc_content=GC_CONTENT)

    # Task 2a and 2b: Generate sequences with errors and track error positions and types
    all_error_data = {}
    for error_rate in ERROR_RATES:
        error_data = generate_sequences_with_error_tracking(template_sequence=dna_sequence,
                                                            error_rate=error_rate,
                                                            num_sequences=NUM_SEQUENCES)
        all_error_data[error_rate] = error_data

    # Task 2c: Plot error position distribution per error type
    plot_error_distribution(all_error_data=all_error_data)

    # Task 3: Plot sequence length distribution per error rate
    plot_length_per_error_rate(all_error_data=all_error_data)

    # Task 4: Plot GC content distribution per error rate
    plot_gc_per_error_rate(all_error_data=all_error_data)

    # Task 6: Plot homopolymer ratios relative to sequence length per nucleotide per error rate
    plot_homopolymer_ratios(all_error_data=all_error_data)

    # Task 6: Plot homopolymer positions per error rate per nucleotide
    plot_homopolymer_positions(all_error_data=all_error_data)


# Execute the main function
if __name__ == '__main__':
    main()
