import requests
import argparse

# NCBI E-Utilities URL
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Function to fetch protein sequences
def fetch_proteins(protein_ids):
    """
    Fetch protein sequences from NCBI using E-Utilities.
    """
    params = {
        "db": "protein",         # Database to query
        "id": ",".join(protein_ids),  # Comma-separated list of protein IDs
        "rettype": "fasta",      # Return type (FASTA format)
        "retmode": "text"        # Return mode (plain text)
    }

    # Make the HTTP request
    response = requests.get(BASE_URL, params=params)

    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching proteins: {response.status_code}")
        return None

# Main function
def main():
    # Argument parser for input and output files
    parser = argparse.ArgumentParser(description="Download protein sequences using NCBI E-Utilities.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing protein IDs (one per line).")
    parser.add_argument("-o", "--output", required=True, help="Output file to save protein sequences.")
    args = parser.parse_args()

    # Read protein IDs from input file
    with open(args.input, "r") as infile:
        protein_ids = [line.strip() for line in infile if line.strip()]

    # Split protein IDs into small batches (recommended by NCBI)
    batch_size = 20
    batches = [protein_ids[i:i + batch_size] for i in range(0, len(protein_ids), batch_size)]

    # Open the output file
    with open(args.output, "w") as outfile:
        for batch in batches:
            print(f"Fetching batch: {', '.join(batch)}")
            sequences = fetch_proteins(batch)
            if sequences:
                outfile.write(sequences + "\n")
            else:
                print(f"Failed to fetch batch: {batch}")

    print(f"Protein sequences saved to {args.output}")

# Run the script
if __name__ == "__main__":
    main()