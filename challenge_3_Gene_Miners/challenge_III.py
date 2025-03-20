import argparse
import pandas as pd
import gzip
import shutil
import os
import urllib.request

def download_goa_file(output_path="goa_human.gaf.gz"):
    '''
    Name:       download_goa_file
    Function:   Downloads the GO annotation file from the EBI FTP server.
    Requires:   No prior download of the file.
    Returns:    A .gz file containing GO annotations.
    '''
    url = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
    
    print(f"Downloading {url}...")
    
    try:
        urllib.request.urlretrieve(url, output_path)
        print("Download complete.")
    except Exception as e:
        print(f"Error downloading the file: {e}")

def unzip_gz_file(gz_path, output_path):
    '''
    Name:       unzip_gz_file
    Function:   Unzips a .gz file to extract its contents.
    Requires:   A valid .gz file path.
    Returns:    A decompressed file ready for use.
    '''
    print(f"Unzipping {gz_path}...")
    with gzip.open(gz_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("Unzipping complete.")

def uniprot_gofinder(uniprot_input, goa_file, output_file):
    '''
    Name:       uniprot_gofinder
    Function:   Retrieves the Top 10 GO terms for each GO branch 
                (Biological Process, Molecular Function, Cellular Component)
                for each UniProt accession number provided in an input file.
    Requires:   - A text file with UniProt accession numbers (one per line).
                - A GO annotation file (downloaded if not available).
    Returns:    A text file with the top 10 GO terms for each category.
    '''
    # Download GOA file if it doesn't exist
    if not os.path.exists(goa_file):
        download_goa_file(goa_file)
    
    # Unzip if necessary
    unzipped_goa_file = goa_file.replace(".gz", "")
    if goa_file.endswith(".gz") and not os.path.exists(unzipped_goa_file):
        unzip_gz_file(goa_file, unzipped_goa_file)

    uniprot_ids = []
    # Read input file to get UniProtIDs
    with open(uniprot_input, 'r') as i_file:
        for line in i_file:
            uniprot_ids.append(line.strip())

    # Read the GOA file into a pandas DataFrame
    df = pd.read_csv(unzipped_goa_file, sep="\t", comment="!", header=None, low_memory=False)

    # Filter rows based on UniProtID
    df = df[df[1].isin(uniprot_ids)]

    # Remove rows where Evidence Code == IEA
    df = df[df[6] != 'IEA']

    # Select UniProtID, GO, and Aspect columns
    df = df.iloc[:, [1, 4, 8]]
    df.columns = ["UniProtID", "GO", "Aspect"]

    # Remove duplicate rows
    df.drop_duplicates(inplace=True)

    # Create dictionaries to store GO terms grouped by category
    bp_terms, mf_terms, cc_terms = {}, {}, {}

    # Write output
    with open(output_file, 'w') as o_file:
        for id in uniprot_ids:
            for _, row in df.iterrows():
                go_term = row["GO"]
                category = row["Aspect"]
                if row["UniProtID"] == id:
                    if category == "P":
                        bp_terms[go_term] = bp_terms.get(go_term, 0) + 1
                    elif category == "F":
                        mf_terms[go_term] = mf_terms.get(go_term, 0) + 1
                    elif category == "C":
                        cc_terms[go_term] = cc_terms.get(go_term, 0) + 1

            # Get the top 10 GO terms by category
            top_10_bp = sorted(bp_terms.items(), key=lambda x: x[1], reverse=True)[:10]
            top_10_mf = sorted(mf_terms.items(), key=lambda x: x[1], reverse=True)[:10]
            top_10_cc = sorted(cc_terms.items(), key=lambda x: x[1], reverse=True)[:10]

            # Get only the GO IDs
            top_10_bp = [go[0] for go in top_10_bp]
            top_10_mf = [go[0] for go in top_10_mf]
            top_10_cc = [go[0] for go in top_10_cc]

            # Write results
            o_file.write(f"UniProtID: {id}\n")
            o_file.write(f"Biological Process: {', '.join(top_10_bp)}\n")
            o_file.write(f"Molecular Function: {', '.join(top_10_mf)}\n")
            o_file.write(f"Cellular Component: {', '.join(top_10_cc)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve the Top 10 GO Terms for each GO branch for each UniProt accession number.")
    parser.add_argument("uniprot_input", type=str, help="Path to file containing a list of UniProt accession numbers (one per line).")
    parser.add_argument("--goa_file", type=str, default="goa_human.gaf.gz", help="Path to GO annotation file (if not provided, it will be downloaded).")
    parser.add_argument("output_file", type=str, help="Path to save the output file.")

    args = parser.parse_args()

    uniprot_gofinder(args.uniprot_input, args.goa_file, args.output_file)
