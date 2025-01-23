import csv

# Define the file paths
file1 = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\PyScriptDelres\\Pyrs333\\mergetables\\tabnew.csv'
file2 = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\PyScriptDelres\\Pyrs333\\mergetables\\tabold.csv'
output_file = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\PyScriptDelres\\Pyrs333\\mergetables\\merged_table.csv'

def read_csv_to_dict(filepath):
    with open(filepath, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        # Strip any leading/trailing whitespace from the headers
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        
        # Debug: print the headers to verify they are correct
        print(f"Headers for {filepath}: {reader.fieldnames}")
        
        data = {row['Sample']: row for row in reader}
    return data, reader.fieldnames

try:
    # Read the first file into a dictionary
    data1, headers1 = read_csv_to_dict(file1)

    # Read the second file into a dictionary
    data2, headers2 = read_csv_to_dict(file2)
except KeyError as e:
    print(f"KeyError: {e}. Make sure the 'Sample' column exists and is correctly spelled in both CSV files.")
    exit(1)

# Get the union of all sample names
all_samples = set(data1.keys()).union(set(data2.keys()))

# Define the headers for the output file
all_headers = ['Sample'] + [header for header in headers1 if header != 'Sample'] + [header for header in headers2 if header != 'Sample']

# Write the merged data to the output file
with open(output_file, mode='w', newline='', encoding='utf-8') as csvfile_out:
    writer = csv.DictWriter(csvfile_out, fieldnames=all_headers, delimiter='\t')
    writer.writeheader()
    
    for sample in all_samples:
        row = {'Sample': sample}
        
        # Check if sample exists in both data1 and data2
        if (sample in data1 and sample in data2 and
            (data1[sample].get('No_filter', '') != data2[sample].get('OLDNo_filter', '') or
             data1[sample].get('Permissive_filter', '') != data2[sample].get('OLDPermissive_filter', '') or
             data1[sample].get('Strict_filter', '') != data2[sample].get('OLDStrict_filter', ''))):
            row['Sample'] += '*'  # Mark with '*'
        
        # Populate data from data1
        if sample in data1:
            row.update({header: data1[sample].get(header, '') for header in headers1 if header != 'Sample'})
        else:
            row.update({header: '' for header in headers1 if header != 'Sample'})
        
        # Populate data from data2
        if sample in data2:
            row.update({header: data2[sample].get(header, '') for header in headers2 if header != 'Sample'})
        else:
            row.update({header: '' for header in headers2 if header != 'Sample'})
        
        writer.writerow(row)

print("Merged file created successfully.")
