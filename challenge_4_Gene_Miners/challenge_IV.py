import argparse

def pathogenic_strains(input_path, output_path):
    '''
    Name:       pathogenic_strains
    Function:   Identifies nucleotide differences between pathogenic and non-pathogenic strains.
    Requires:   A FASTA file containing sequences labeled as "pathogenic" or "non_pathogenic".
    Returns:    A text file listing nucleotide positions where differences occur between groups.
    '''
    # Read input file to get FASTA sequences
    with open(input_path, 'r') as i_file:
        sequences = {}
        while True:
            # We read each line
            line = i_file.readline()
            # eliminate \n
            line = line.replace("\n", "")
            # exit if there are no more lines left
            if not line:
                break
            # if it starts by '>', we are in front of a new sequence
            elif line[0] == '>':
                # saving the label
                sequence_label = line[1:]
            # if not, it is a sequence
            else:
                # update dictionary using {label:sequence} as {key:value}
                sequences.update({sequence_label:line})

        # Now we need to see wich nucleotides are common between the strains in the 2 groups
        pathogenic = {}
        non_pathogenic = {}

        # classify each key by his name
        for key in sequences.keys():
            if key.startswith("pathogenic"):
                pathogenic.update({key:sequences[key]})
            elif key.startswith("non_pathogenic"):
                non_pathogenic.update({key:sequences[key]})


        ### compare each group's sequences
        values = list(pathogenic.values())

        # Write output
        with open(output_path, 'w') as o_file:
            # for each nucleotide
            for i in range(len(values[0])):
                # reset sets
                pat_variants = set()
                nonpat_variants = set()
                # for each sequence in pathogenic group
                for sequence in pathogenic.values():
                    # save the variant for the specific nucleotide
                    pat_variants.add(sequence[i])
                # for each sequence in non-pathogenic group
                for sequence in non_pathogenic.values():
                    # save the variant for the specific nucleotide
                    nonpat_variants.add(sequence[i])
                
                # if the sets don't share any element (are disjoint)
                if pat_variants.isdisjoint(nonpat_variants) == True:
                    
                    # return the differences
                    o_file.write("Position {}: ".format(i+1) + " Pathogenic -> {}".format(pat_variants) + ", Non-Pathogenic -> {}".format(nonpat_variants) + "\n")

if __name__ == "__main__":
    # Argument parser setup for command-line execution
    parser = argparse.ArgumentParser(description="Identifies nucleotide differences between pathogenic and non-pathogenic strains.")
    parser.add_argument("input_file", type=str, help="Path to input file")
    parser.add_argument("output_file", type=str, help="Path to output file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Execute the main function
    pathogenic_strains(args.input_file, args.output_file)
