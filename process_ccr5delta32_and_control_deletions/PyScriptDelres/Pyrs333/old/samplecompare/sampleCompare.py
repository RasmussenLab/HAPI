import pandas as pd

# File paths for input files
newfile_path = r'C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\PyScriptDelres\Pyrs333\samplecompare\New.tsv'
oldfile_path = r'C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\PyScriptDelres\Pyrs333\samplecompare\Old.tsv'

# Read newfile
newfile = pd.read_csv(newfile_path, delimiter='\t')

# Read oldfile
oldfile = pd.read_csv(oldfile_path, delimiter='\t')

# Merge DataFrames on Assigned and Sample columns with an outer join to include all samples
merged = pd.merge(newfile, oldfile, on=['Assigned', 'Sample'], suffixes=('_new', '_old'), how='outer')

# Function to mark differences
def mark_difference(row):
    if pd.isnull(row['No_filter_new']) and pd.notnull(row['No_filter_old']):
        return '*'  # Mark if sample is only in oldfile
    elif pd.isnull(row['No_filter_old']) and pd.notnull(row['No_filter_new']):
        return '*'  # Mark if sample is only in newfile
    elif (row['No_filter_new'] != row['No_filter_old'] or
          row['Permissive_filter_new'] != row['Permissive_filter_old'] or
          row['Strict_filter_new'] != row['Strict_filter_old']):
        return '*'  # Mark if there are differences
    else:
        return ''   # No differences or missing data

# Add Difference column
merged['Difference'] = merged.apply(mark_difference, axis=1)

# Drop rows with all NaN values (if any)
merged.dropna(how='all', inplace=True)

# Output file path
output_file_path = r'C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\PyScriptDelres\Pyrs333\samplecompare\output.tsv'

# Save the merged DataFrame to a new file
merged.to_csv(output_file_path, sep='\t', index=False)

print(f"Output saved to {output_file_path}")
