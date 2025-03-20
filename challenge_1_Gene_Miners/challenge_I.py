import argparse
import numpy as np

def import_value(line, parameter_name):
    '''
    Name:       import_value
    Function:   gets the value of a parameter found in a line of a text.
    Requires:   a string ("line") starting with the name of the parameter,
                the "parameter_name" to find in the string.
    Returns:    a string containing the value of the parameter.
    '''
    # eliminate spaces
    line = line.replace(" ", "")
    # eliminate prefix
    line = line.removeprefix(parameter_name + "=")

    return line

def mutation_rate_bootstrap(input_file, output_file):
    '''
    Name:       mutation_rate_bootstrap
    Function:   Calculates mutation rate per site per generation (MR)
                and the corresponding 95% confidence interval.
    Requires:   A .txt file containing:
                - A vector "mut" of comma-separated integers representing
                  the total number of mutations observed in each sampled individual.
                - An integer "l" representing the genome length in base pairs.
                - An integer "g" representing the total number of generations.
                - An integer "b" representing the number of bootstrap replicates requested.
    Returns:    A .txt file containing:
                - The mean mutation rate (MR).
                - A 95% confidence interval calculated using bootstrapping.
    '''

    # Initialize variables
    mut = []
    l = 0
    g = 0
    b = 0

    # Read input file
    with open(input_file, 'r') as i_f:
        for line in i_f:
            if line.startswith("mut"):
                # Extract mutation vector and convert to integer list
                mut_line = import_value(line, "mut")
                mut = list(map(int, mut_line.split(",")))
            elif line.startswith("l"):
                # Extract genome length
                l = int(import_value(line, "l"))
            elif line.startswith("g"):
                # Extract total number of generations
                g = int(import_value(line, "g"))
            elif line.startswith("b"):
                # Extract number of bootstrap replicates
                b = int(import_value(line, "b"))

    # Compute the observed mutation rate (MR)
    mean_mutation_rate = sum(mut) / (l * g)

    # Bootstrap resampling to estimate confidence interval
    bootstrap_rates = []
    for _ in range(b):
        # Sample with replacement from the original mutation data
        sample = np.random.choice(mut, size=len(mut), replace=True)
        bootstrap_rate = sum(sample) / (l * g)
        bootstrap_rates.append(bootstrap_rate)

    # Calculate 95% confidence interval (2.5th and 97.5th percentiles)
    lower_bound = np.percentile(bootstrap_rates, 2.5)
    upper_bound = np.percentile(bootstrap_rates, 97.5)

    # Write output file
    with open(output_file, "w") as o_f:
        o_f.write(f"Mean mutation rate observed: {mean_mutation_rate:.2e}\n")
        o_f.write(f"95% Confidence Interval: [{lower_bound:.1e}, {upper_bound:.1e}]\n")

if __name__ == "__main__":
    # Argument parser setup for command-line execution
    parser = argparse.ArgumentParser(description="Calculates mutation rate per site per generation (MR) and the corresponding 95% confidence interval using bootstrap resampling.")
    parser.add_argument("input_file", type=str, help="Path to input file")
    parser.add_argument("output_file", type=str, help="Path to output file")

    # Parse command-line arguments
    args = parser.parse_args()

    # Execute the main function
    mutation_rate_bootstrap(args.input_file, args.output_file)
