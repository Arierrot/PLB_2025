import argparse

def gc_content(sequence):
    '''
    Name:       gc_content
    Function:   Calculates the GC content percentage in a given DNA sequence.
    Requires:   A string representing a DNA sequence.
    Returns:    A float representing the GC content percentage.
    '''
    gc_count = sequence.count('G') + sequence.count('C')
    return round((gc_count / len(sequence)) * 100, 1)

def longest_common_substring(seq1, seq2):
    '''
    Name:       longest_common_substring
    Function:   Finds the longest common substring between two sequences using dynamic programming.
    Requires:   Two strings representing DNA sequences.
    Returns:    A string representing the longest common substring.
    '''
    m = len(seq1) 
    n = len(seq2)

    # Creating a 2D matrix with 0s (size = len(seq1)+1 x len(seq2)+1)
    matrix = [[0] * (n+1) for _ in range(m+1)]

    max_len = 0 # Saving the length of the Longest Common Substring (LCS)
    end_pos = 0 # Saving the ending position where LCS ends in seq1

    # Filling the matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:  # If the characters match...
                matrix[i][j] = matrix[i-1][j-1] + 1  # Add 1 to the diagonal
    
                # Refreshing the max subsequence found
                if matrix[i][j] > max_len:
                    max_len = matrix[i][j]
                    end_pos = i  # Ending position where LCS ends in seq1
    
    # Return Longest Common Substring
    lcs = seq1[end_pos - max_len : end_pos]

    return lcs

def process_sequences(input_file, output_file):
    '''
    Name:       process_sequences
    Function:   Reads an input file containing two DNA sequences, calculates GC content,
                sequence lengths, and finds the longest common substring.
    Requires:   - A text file with two DNA sequences.
    Returns:    A text file with calculated sequence lengths, GC content percentages,
                and the longest common substring.
    '''
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    seq1 = lines[0].strip().split()[1]
    seq2 = lines[1].strip().split()[1]

    # Compute values
    len_seq1, len_seq2 = len(seq1), len(seq2)
    gc_seq1, gc_seq2 = gc_content(seq1), gc_content(seq2)
    lcs = longest_common_substring(seq1, seq2)

    # Write output file
    with open(output_file, 'w') as f:
        f.write(f"seq1 length: {len_seq1}\n")
        f.write(f"seq2 length: {len_seq2}\n")
        f.write(f"seq1 GC content: {gc_seq1}%\n")
        f.write(f"seq2 GC content: {gc_seq2}%\n")
        f.write(f"Longest common substring: {lcs}\n")

if __name__ == "__main__":
    '''
    Name:       Main Execution
    Function:   Handles command-line arguments and calls process_sequences function.
    Requires:   - Input and output file paths as command-line arguments.
    Returns:    Generates an output file with sequence analysis results.
    '''
    # Argument parser setup for command-line execution
    parser = argparse.ArgumentParser(description="Analyze aptamer sequences for GC content and longest common substring.")
    parser.add_argument("input_file", type=str, help="Path to input file")
    parser.add_argument("output_file", type=str, help="Path to output file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Execute the main function
    process_sequences(args.input_file, args.output_file)
