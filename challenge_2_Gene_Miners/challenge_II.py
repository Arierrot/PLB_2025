import argparse
from Bio import Entrez

def search_gene_names(organism, output_file):
    '''
    Name:       search_gene_names
    Function:   Queries NCBI Entrez API to retrieve a list of gene names
                for the specified organism and saves them to a file.
    Requires:   - organism (str): Scientific name of the organism.
                - output_file (str): Path to the output file to store gene names.
    Returns:    A .txt file containing gene names (one per line).
    '''

    # Set Entrez email (Required by NCBI), just for notifications
    Entrez.email = "anonymous@gmail.com"

    # Search for the organism's genes in the NCBI Gene database
    search_term = f"{organism}[Organism]"
    handle = Entrez.esearch(db="gene", term=search_term)
    record = Entrez.read(handle)
    handle.close()

    # Extract gene IDs
    gene_ids = record["IdList"]

    # Fetch detailed gene information
    gene_names = []
    for gene_id in gene_ids:
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        record = Entrez.read(handle)
        handle.close()

        try:
            # Extract official gene symbol
            gene_symbol = record[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
            gene_names.append(gene_symbol)
        except KeyError:
            continue  # Skip if there's no gene symbol

    # Write gene names to output file
    with open(output_file, "w") as f:
        for gene in gene_names:
            f.write(f"{gene}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve gene list from NCBI Entrez API.")
    parser.add_argument("organism", type=str, help="Scientific name of the organism (e.g., 'Homo sapiens').")
    parser.add_argument("output_file", type=str, help="Path to save the output file with gene names.")

    # Parse command-line arguments
    args = parser.parse_args()

    # Execute the main function
    search_gene_names(args.organism, args.output_file)