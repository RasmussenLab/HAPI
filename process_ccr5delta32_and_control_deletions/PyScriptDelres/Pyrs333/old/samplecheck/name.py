import os
import pandas as pd

# Get the current script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Specify the path to input.tsv relative to the script's location
input_file = os.path.join(script_dir, "input.tsv")

# Check if the file exists
if not os.path.exists(input_file):
    print(f"Error: File '{input_file}' not found. Please check the file path.")
    exit(1)

# Read the input TSV file
try:
    df = pd.read_csv(input_file, sep='\t')

    # Extract sample names from FTP URLs
    df['Sample Name'] = df['FTP_URL'].apply(lambda url: url.split("/")[-1].replace(".final.bam", ""))

    # Write sample names to mylist.tsv
    output_file = os.path.join(script_dir, 'mylist.tsv')
    df[['Sample Name']].to_csv(output_file, sep='\t', index=False)

    print(f"Extracted sample names saved to {output_file}.")
except Exception as e:
    print(f"Error occurred while processing file: {e}")
