import csv

# Define the file paths
file1 = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\New_files150724\\PyScriptDelres\\Pyrs333\\mergetables\\tabnew.csv'
file2 = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\New_files150724\\PyScriptDelres\\Pyrs333\\mergetables\\tabold.csv'
output_file = 'C:\\Users\\k_rav\\OneDrive\\Skrivebord\\Speciale\\1NOV22120\\Cell\\Cell_reviewers_comments\\Deletions_background\\Results_background\\New_files150724\\PyScriptDelres\\Pyrs333\\mergetables\\merged_table.csv'

def read_csv_to_dict(filepath):
    with open(filepath, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        print(f"Headers for {filepath}: {reader.fieldnames}")
        data = {row['Sample']: row for row in reader}
    return data, reader.fieldnames

try:
    data1, headers1 = read_csv_to_dict(file1)
    data2, headers2 = read_csv_to_dict(file2)
except KeyError as e:
    print(f"KeyError: {e}. Make sure the 'Sample' column exists and is correctly spelled in both CSV files.")
    exit(1)

all_samples = set(data1.keys()).union(set(data2.keys()))
all_headers = ['Sample'] + [header for header in headers1 if header != 'Sample'] + [header for header in headers2 if header != 'Sample']

starred_samples = []
blank_samples = []
regular_samples = []

for sample in all_samples:
    row = {'Sample': sample}
    
    if sample in data1:
        row.update({header: data1[sample].get(header, '') for header in headers1 if header != 'Sample'})
    else:
        row.update({header: '' for header in headers1 if header != 'Sample'})
    
    if sample in data2:
        row.update({header: data2[sample].get(header, '') for header in headers2 if header != 'Sample'})
    else:
        row.update({header: '' for header in headers2 if header != 'Sample'})
    
    if (sample in data1 and sample in data2 and
        (data1[sample].get('No_filter', '') != data2[sample].get('OLDNo_filter', '') or
         data1[sample].get('Permissive_filter', '') != data2[sample].get('OLDPermissive_filter', '') or
         data1[sample].get('Strict_filter', '') != data2[sample].get('OLDStrict_filter', ''))):
        row['Sample'] += '*'
        starred_samples.append(row)
    elif '' in row.values():
        blank_samples.append(row)
    else:
        regular_samples.append(row)

with open(output_file, mode='w', newline='', encoding='utf-8') as csvfile_out:
    writer = csv.DictWriter(csvfile_out, fieldnames=all_headers, delimiter='\t')
    writer.writeheader()
    
    for row in starred_samples:
        writer.writerow(row)
    
    for row in blank_samples:
        writer.writerow(row)
    
    for row in regular_samples:
        writer.writerow(row)

print("Merged file created successfully.")


